"""
Génération d'un fichier résultat 2D (format Telemac) avec plusieurs variables et éventuellement plusieurs temps.
La taille des éléments du maillage est controlable.
Seulement les résultats (ou données) aux sections rattachées à des branches 20 sont considérées.

Si un fichier rcal est spécifié alors on peut traiter :
- tous les calculs permanents
- un calcul transitoire en spécifiant son nom dans l'argument `--calc_trans`

Les variables écrites dans le fichier de sortie sont :
* FOND
* IS LIT ACTIVE (0 = lit inactif, 1 = lit actif)
* FROTTEMENT (moyenne sur la verticale si plusieurs valeurs)
* Les variables supplémentaires si un résultat est lu sont :
    * HAUTEUR D'EAU (la variable 'Z' est nécessaire)
    * VITESSE SCALAIRE (seulement si la variable 'Vact' est présente)

TODO:
* Tenir compte des lignes de contrainte (lit numérotés)
* Support des casiers?
"""
import numpy as np
import os.path
from pyteltools.slf import Serafin
import sys
from time import perf_counter

from crue10.emh.branche import Branche
from crue10.emh.section import SectionProfil
from crue10.run import CrueRun
from crue10.study import Study
from crue10.utils import CrueError

from core.arg_command_line import MyArgParse
from core.base import LigneContrainte, MeshConstructor, ProfilTravers, SuiteProfilsTravers
from core.utils import logger, set_logger_level, TatooineException


LANG = 'fr'  # for variable names
VARIABLES_FROM_GEOMETRY = ['FOND', 'IS LIT ACTIVE', 'FROTTEMENT']


def mesh_crue10_run(args):
    set_logger_level(args.verbose)
    t1 = perf_counter()

    # Read the model and its submodels from xml/shp files
    study = Study(args.infile_etu)
    model = study.get_model(args.model_name)
    model.read_all()
    logger.info(model)
    for submodel in model.submodels:
        submodel.remove_sectioninterpolee()
        submodel.normalize_geometry()
        logger.info(submodel.summary())
        # submodel.write_shp_limites_lits_numerotes('limites_lits.shp')  # DEBUG
    logger.info(model)

    triangles = {}
    mesh_constr = None
    points = None
    suite = []

    id_profile = 0
    for i, branche in enumerate(model.get_branche_list()):
        # Ignore branch if branch_patterns is set and do not match with current branch name
        if args.branch_patterns is not None:
            ignore = True
            for pattern in args.branch_patterns:
                if pattern in branche.id:
                    ignore = False
                    break
        else:
            ignore = False

        if branche.type != 20:
            ignore = True

        if not ignore:
            logger.info("===== TRAITEMENT DE LA BRANCHE %s =====" % branche.id)
            axe = branche.geom
            try:
                profils_travers = SuiteProfilsTravers()
                for section in branche.sections:
                    if isinstance(section, SectionProfil):
                        coords = list(section.get_coord(add_z=True))
                        profile = ProfilTravers(section.id, [(coord[0], coord[1]) for coord in coords], 'Section')

                        # Determine some variables (constant over the simulation) from the geometry
                        z = np.array([coord[2] for coord in coords])
                        is_bed_active = section.get_is_bed_active_array()
                        mean_strickler = section.get_friction_coeff_array()
                        profile.coord.values = np.core.records.fromarrays(
                            np.column_stack((z, is_bed_active, mean_strickler)).T,
                            names=VARIABLES_FROM_GEOMETRY
                        )

                        profils_travers.add_profile(profile)
                        id_profile += 1

                profils_travers.compute_dist_proj_axe(axe, args.dist_max)
                if len(profils_travers) >= 2:
                    profils_travers.check_intersections()
                    # profils_travers.sort_by_dist() is useless because profiles are already sorted
                    lignes_contraintes = LigneContrainte.get_lines_and_set_limits_from_profils(profils_travers)

                    mesh_constr = MeshConstructor(profils_travers, args.pas_trans, args.nb_pts_trans)
                    mesh_constr.build_interp(lignes_contraintes, args.pas_long, args.constant_ech_long)

                    mesh_constr.build_mesh(True)
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
                    logger.warning("Branche ignorée par manque de sections")
            except TatooineException as e:
                logger.error("/!\ Branche ignorée à cause d'une erreur bloquante :")
                logger.error(e.message)
            logger.info("\n")

    if points is None:
        raise CrueError("Aucun point à traiter, adaptez l'option `--branch_patterns`")

    mesh_constr.points = points
    mesh_constr.profils_travers = suite
    mesh_constr.triangle = triangles

    logger.info(mesh_constr.summary())  # General information about the merged mesh

    if args.infile_rcal:
        # Read rcal result file
        run = CrueRun(args.infile_rcal)
        logger.info(run.summary())

        # Check result consistency
        missing_sections = model.get_missing_active_sections(run.emh['Section'])
        if missing_sections:
            raise CrueError("Sections manquantes :\n%s" % missing_sections)

        # Subset results to get requested variables at active sections
        varnames_1d = run.variables['Section']
        logger.info("Variables 1D disponibles aux sections: %s" % varnames_1d)
        pos_z = varnames_1d.index('Z')
        pos_variables = [run.variables['Section'].index(var) for var in varnames_1d]
        pos_sections_list = [run.emh['Section'].index(profil.id) for profil in mesh_constr.profils_travers]

        additional_variables_id = ['H']
        if 'Vact' in varnames_1d:
            additional_variables_id.append('M')

        values_geom = mesh_constr.interp_values_from_profiles()
        z_bottom = values_geom[0, :]
        with Serafin.Write(args.outfile_mesh, LANG, overwrite=True) as resout:
            output_header = Serafin.SerafinHeader(title='%s (written by tatooinemesher)' % os.path.basename(args.outfile_mesh),
                                                  lang=LANG)
            output_header.from_triangulation(mesh_constr.triangle['vertices'],
                                             mesh_constr.triangle['triangles'] + 1)
            for var_name in VARIABLES_FROM_GEOMETRY:
                output_header.add_variable_str(var_name, var_name, '')
            for var_id in additional_variables_id:
                output_header.add_variable_from_ID(var_id)
            for var_name in varnames_1d:
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
                    values_1d = mesh_constr.interp_from_values_at_profiles(variables_at_profiles)

                    # Compute water depth: H = Z - Zf and clip below 0m (avoid negative values)
                    depth = np.clip(values_1d[pos_z, :] - z_bottom, a_min=0.0, a_max=None)

                    # Merge and write values
                    if 'Vact' in varnames_1d:
                        # Compute velocity magnitude from Vact and apply mask "is active bed"
                        velocity = values_1d[varnames_1d.index('Vact'), :] * values_geom[1, :]
                        values = np.vstack((values_geom, depth, velocity, values_1d))
                    else:
                        values = np.vstack((values_geom, depth, values_1d))
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
                    values_1d = mesh_constr.interp_from_values_at_profiles(variables_at_profiles)

                    # Compute water depth: H = Z - Zf and clip below 0m (avoid negative values)
                    depth = np.clip(values_1d[pos_z, :] - z_bottom, a_min=0.0, a_max=None)

                    # Merge and write values
                    if 'Vact' in varnames_1d:
                        # Compute velocity magnitude from Vact and apply mask "is active bed"
                        velocity = values_1d[varnames_1d.index('Vact'), :] * values_geom[1, :]
                        values = np.vstack((values_geom, depth, velocity, values_1d))
                    else:
                        values = np.vstack((values_geom, depth, values_1d))
                    resout.write_entire_frame(output_header, time, values)

    else:
        # Write a single frame with only variables from geometry
        mesh_constr.export_mesh(args.outfile_mesh)

    t2 = perf_counter()
    logger.info("=> le temps d'execution est de : {}s".format(t2-t1))


parser = MyArgParse(description=__doc__)
# Inputs
parser_infiles = parser.add_argument_group("Modèle et run Crue10 à traiter")
parser_infiles.add_argument("infile_etu", help="fichier d'étude Crue10 (etu.xml)")
parser_infiles.add_argument("model_name", help="nom du modèle")
parser_infiles.add_argument("--infile_rcal", help="fichier de résultat (rcal.xml)")
parser_infiles.add_argument("--calc_trans", help="nom du calcul transitoire à traiter "
                                                 "(sinon considère tous les calculs permanents)")
# Branch parameters
parser_branch = parser.add_argument_group("Paramètres pour choisir les branches à traiter")
parser_branch.add_argument("--branch_types_filter", nargs='+', default=Branche.TYPES_WITH_LENGTH,
                            help="types des branches à traiter")
parser_branch.add_argument("--branch_patterns", nargs='+', default=None,
                            help="chaîne(s) de caractères pour ne conserver que les branches dont le nom contient"
                                 " une de ces expressions")
# Mesh parameters
parser_mesh = parser.add_argument_group("Paramètres pour la génération du maillage 2D")
parser_mesh.add_argument("--dist_max", type=float, help="distance de recherche maxi des 'intersections fictifs' "
                                                        "pour les limites de lits (en m)", default=0.01)
parser.add_argument("--pas_long", type=float, help="pas d'interpolation longitudinal (en m)")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--pas_trans", type=float, help="pas d'interpolation transversal (en m)")
group.add_argument("--nb_pts_trans", type=int, help="nombre de noeuds transveralement")
parser_mesh.add_argument("--constant_ech_long",
                         help="méthode de calcul du nombre de profils interpolés entre profils : "
                              "par profil (constant, ie True) ou par lit (variable, ie False)",
                         action='store_true')
# Outputs
parser_outfiles = parser.add_argument_group('Fichier de sortie')
parser_outfiles.add_argument("outfile_mesh", help="Résultat 2D au format slf (Telemac)")

if __name__ == '__main__':
    args = parser.parse_args()
    try:
        mesh_crue10_run(args)
    except FileNotFoundError as e:
        logger.critical(e)
        sys.exit(1)
    except CrueError as e:
        logger.critical(e)
        sys.exit(2)
