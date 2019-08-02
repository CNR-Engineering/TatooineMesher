"""
Génération d'un résultat 2D (format Telemac) avec plusieurs variables et plusieurs temps
La taille des éléments du maillage est controlable.
Only results at sections (type 20) are considered.
Consider either all steady state calculations or a single transient calculation (specify `--calc_trans` argument).

TODO:
* Add constrainted lines between "lit numérotés"
* Add Vact with filter on active beds
* Compute Strickler coefficient (called `Kmin` or `Kact_eq` in rcal) and TAU (called `Tauact`)
* Handle casier?
"""
import numpy as np
import os.path
from pyteltools.slf import Serafin
import sys
from time import perf_counter

from crue10.emh.section import SectionProfil
from crue10.run import CrueRun
from crue10.study import Study
from crue10.utils import CrueError

from core.arg_command_line import MyArgParse
from core.base import LigneContrainte, MeshConstructor, ProfilTravers, SuiteProfilsTravers
from core.utils import float_vars, logger, set_logger_level, TatooineException


LANG = 'fr'  # for variable names
VARIABLES = ['Z', 'Fr']  # First variable has to be Z
ADDITIONAL_VARIABLES = ['Zf', 'H']  # Hardcoded below


def mesh_crue10_run(args):
    set_logger_level(args.verbose)
    t1 = perf_counter()

    # Read the model and its submodels from xml/shp files
    try:
        study = Study(args.infile_etu)
        model = study.get_model(args.model_name)
        model.read_all()
        logger.info(model)
        for submodel in model.submodels:
            submodel.remove_sectioninterpolee()
            submodel.normalize_geometry()
            logger.info(submodel.summary())
    except FileNotFoundError as e:
        logger.critical(e)
        sys.exit(1)
    except CrueError as e:
        logger.critical(e)
        sys.exit(2)
    logger.info(model)

    triangles = {}
    mesh_constr = None
    points = None
    suite = []

    id_profile = 0
    for i, branche in enumerate(model.get_branche_list()):
        if branche.type == 20:
            logger.info("===== TRAITEMENT DE LA BRANCHE %s =====" % branche.id)
            axe = branche.geom
            try:
                profils_travers = SuiteProfilsTravers()
                for section in branche.sections:
                    if isinstance(section, SectionProfil):
                        coords = list(section.get_coord(add_z=True))
                        profile = ProfilTravers(section.id, [(coord[0], coord[1]) for coord in coords], 'Section')
                        profile.coord.values = np.array([(coord[2],) for coord in coords], dtype=float_vars('Z'))
                        profils_travers.add_profile(profile)
                        id_profile += 1

                if len(profils_travers) >= 2:
                    profils_travers.compute_dist_proj_axe(axe, args.dist_max)
                    profils_travers.check_intersections()
                    # profils_travers.sort_by_dist() is useless because profiles are already sorted
                    lignes_contraintes = LigneContrainte.get_lines_from_profils(profils_travers)
                    profils_travers.find_and_add_limits(lignes_contraintes, args.dist_max)

                    mesh_constr = MeshConstructor(profils_travers, args.pas_trans)
                    mesh_constr.build_interp(lignes_contraintes, args.pas_long, args.constant_ech_long)

                    mesh_constr.build_mesh()
                    suite += mesh_constr.profils_travers
                    if not triangles:  # set initial values from first iteration
                        points = mesh_constr.points
                        triangles['triangles'] = mesh_constr.triangle['triangles']
                        triangles['vertices'] = mesh_constr.triangle['vertices']
                    else:  # concatenate with current sub-mesh for next iterations
                        last_zone = points['zone'].max()
                        mesh_constr.points['zone'] += last_zone + 2
                        points = np.hstack((points, mesh_constr.points))
                        triangles['triangles'] = np.vstack(
                            (triangles['triangles'],
                             triangles['vertices'].shape[0] + mesh_constr.triangle['triangles']))
                        triangles['vertices'] = np.vstack((triangles['vertices'], mesh_constr.triangle['vertices']))
                else:
                    logger.info("Branche ignorée par manque de sections")
            except TatooineException as e:
                logger.critical("/!\ Branche ignorée à cause d'une erreur bloquante :")
                logger.critical(e.message)
            logger.info("\n")

    mesh_constr.points = points
    mesh_constr.profils_travers = suite
    mesh_constr.triangle = triangles

    logger.info(mesh_constr.summary())  # General information about the merged mesh

    # Read rcal result file
    run = CrueRun(args.infile_rcal)
    logger.info(run.summary())

    # Check result consistency
    missing_sections = model.get_missing_active_sections(run.emh['Section'])
    if missing_sections:
        logger.error("Missing sections:\n%s" % missing_sections)
        return

    # Subset results to get requested variables at active sections
    pos_variables = [run.variables['Section'].index(var) for var in VARIABLES]
    pos_sections_list = [run.emh['Section'].index(profil.id) for profil in mesh_constr.profils_travers]

    z_bottom = mesh_constr.interp_values_from_profiles()[0, :]
    with Serafin.Write(args.outfile_mesh, LANG, overwrite=True) as resout:
        output_header = Serafin.SerafinHeader(title='%s (written by tatooinemesher)' % os.path.basename(args.outfile_mesh),
                                              lang=LANG)
        output_header.from_triangulation(mesh_constr.triangle['vertices'],
                                         mesh_constr.triangle['triangles'] + 1)
        for var_name in VARIABLES + ADDITIONAL_VARIABLES:
            output_header.add_variable_str(var_name, var_name, '')
        resout.write_header(output_header)

        if args.calc_trans is None:
            for i, calc_name in enumerate(run.calc_perms.keys()):
                logger.info("~> Calcul permanent %s" % calc_name)
                # Read a single *steady* calculation
                res_perm = run.get_res_perm(calc_name)
                res = res_perm['Section']
                variables_at_profiles = res[pos_sections_list, :][:, pos_variables]

                # Interpolate between sections
                values = mesh_constr.interp_from_values_at_profiles(variables_at_profiles)

                # Add additional variables
                # H = water depth = Z - Zf and clip below 0m (avoid negative values)
                depth = np.clip(values[0, :] - z_bottom, a_min=0.0, a_max=None)

                # Merge and write values
                values = np.vstack((values, z_bottom, depth))
                resout.write_entire_frame(output_header, 3600.0 * i, values)
        else:
            res_trans = run.get_calc_trans(args.calc_trans)
            logger.info("Calcul transitoire %s" % args.calc_trans)
            res_all = run.get_res_trans(args.calc_trans)['Section']

            for i, (time, _) in enumerate(res_trans.frame_list):
                logger.info("~> %fs" % time)
                res = res_all[i, :, :]
                variables_at_profiles = res[pos_sections_list, :][:, pos_variables]

                # Interpolate between sections
                values = mesh_constr.interp_from_values_at_profiles(variables_at_profiles)

                # Add additional variables
                # H = water depth = Z - Zf and clip below 0m (avoid negative values)
                depth = np.clip(values[0, :] - z_bottom, a_min=0.0, a_max=None)

                # Merge and write values
                values = np.vstack((values, z_bottom, depth))
                resout.write_entire_frame(output_header, time, values)

    t2 = perf_counter()
    logger.info("=> le temps d'execution est de : {}s".format(t2-t1))


parser = MyArgParse(description=__doc__)
# Inputs
parser_infiles = parser.add_argument_group("Modèle et run Crue10 à traiter")
parser_infiles.add_argument("infile_etu", help="fichier d'étude Crue10 (etu.xml)")
parser_infiles.add_argument("model_name", help="nom du modèle")
parser_infiles.add_argument("infile_rcal", help="fichier de résultat (rcal.xml)")
parser_infiles.add_argument("--calc_trans", help="nom du calcul transitoire à traiter "
                                                 "(sinon considère tous les calculs permanents)")
# Mesh parameters
parser_mesh = parser.add_argument_group("Paramètres pour la génération du maillage 2D")
parser_mesh.add_argument("--pas_long", type=float, help="pas d'interpolation longitudinal (en m)", default=5)
parser_mesh.add_argument("--pas_trans", type=float, help="pas d'interpolation transversal (en m)", default=3.5)
parser_mesh.add_argument("--dist_max", type=float, help="distance de recherche maxi des 'intersections fictifs' "
                                                        "pour les limites de lits (en m)", default=0.01)
parser_mesh.add_argument("--constant_ech_long",
                         help="méthode de calcul du nombre de profils interpolés entre profils : "
                              "par profil (constant, ie True) ou par lit (variable, ie False)",
                         action='store_true')
# Outputs
parser_outfiles = parser.add_argument_group('Fichier de sortie')
parser_outfiles.add_argument("--outfile_mesh", help="Résultat 2D au format slf (Telemac)")

if __name__ == '__main__':
    args = parser.parse_args()
    mesh_crue10_run(args)
