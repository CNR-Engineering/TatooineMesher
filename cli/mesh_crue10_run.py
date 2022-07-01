"""
mesh_crue10_run.py

Génération d'un fichier résultat 2D (format Telemac) avec plusieurs variables et éventuellement plusieurs temps.
La taille des éléments du maillage est controlable.

Si un fichier rcal est spécifié alors on peut traiter :
- tous les calculs permanents
- un calcul transitoire en spécifiant son nom dans l'argument `--calc_unsteady`

Les variables écrites dans le fichier de sortie sont :
* FOND
* IS LIT ACTIVE (0 = lit inactif, 1 = lit actif)
* FROTTEMENT (moyenne sur la verticale si plusieurs valeurs)
* Les variables supplémentaires si un résultat est lu sont :
    * HAUTEUR D'EAU (la variable 'Z' aux sections/casiers est nécessaire.
        Attention, il ne faut pas avoir sorti la charge 'H' aux sections)
    * VITESSE SCALAIRE (seulement si la variable 'Vact' est présente)

Seulement les branches et les casiers actifs sont traités.
"""
import numpy as np
from numpy.lib.recfunctions import unstructured_to_structured
import os.path
from pyteltools.slf import Serafin
import sys
from time import perf_counter
import triangle

from crue10.emh.branche import Branche
from crue10.emh.section import SectionProfil
from crue10.etude import Etude
from crue10.run.results import RunResults
from crue10.utils import ExceptionCrue10

from tatooinemesher.constraint_line import ConstraintLine
from tatooinemesher.mesh_constructor import MeshConstructor
from tatooinemesher.section import CrossSection, CrossSectionSequence
from tatooinemesher.interp.raster import interp_raster
from tatooinemesher.utils.arg_command_line import MyArgParse
from tatooinemesher.utils import logger, resample_2d_line, set_logger_level, TatooineException


VARIABLES_FROM_GEOMETRY = ['B', 'IS LIT ACTIVE', 'W']


def mesh_crue10_run(args):
    set_logger_level(args.verbose)
    t1 = perf_counter()

    # Read the model and its submodels from xml/shp files
    etude = Etude(args.infile_etu)
    modele = etude.get_modele(args.model_name)
    modele.read_all()
    logger.info(modele)
    for sous_modele in modele.liste_sous_modeles:
        sous_modele.remove_sectioninterpolee()
        sous_modele.normalize_geometry()
        logger.info(sous_modele.summary())
        # sous_modele.write_shp_limites_lits_numerotes('limites_lits.shp')  # DEBUG
    logger.info(modele)

    global_mesh_constr = MeshConstructor()

    # Handle branches in minor bed
    for i, branche in enumerate(modele.get_liste_branches()):
        # Ignore branch if branch_patterns is set and do not match with current branch name
        if args.branch_patterns is not None:
            ignore = True
            for pattern in args.branch_patterns:
                if pattern in branche.id:
                    ignore = False
                    break
        else:
            ignore = False

        if branche.type not in args.branch_types_filter or not branche.is_active:
            ignore = True

        if not ignore:
            logger.info("===== TRAITEMENT DE LA BRANCHE %s =====" % branche.id)
            axe = branche.geom
            try:
                section_seq = CrossSectionSequence()
                for crue_section in branche.liste_sections_dans_branche:
                    if isinstance(crue_section, SectionProfil):
                        coords = list(crue_section.get_coord(add_z=True))
                        section = CrossSection(crue_section.id, [(coord[0], coord[1]) for coord in coords], 'Section')

                        # Determine some variables (constant over the simulation) from the geometry
                        z = np.array([coord[2] for coord in coords])
                        is_bed_active = crue_section.get_is_bed_active_array()
                        mean_strickler = crue_section.get_friction_coeff_array()
                        section.coord.values = np.core.records.fromarrays(
                            np.column_stack((z, is_bed_active, mean_strickler)).T,
                            names=VARIABLES_FROM_GEOMETRY
                        )

                        section_seq.add_section(section)

                section_seq.compute_dist_proj_axe(axe, args.dist_max)
                if len(section_seq) >= 2:
                    section_seq.check_intersections()
                    # section_seq.sort_by_dist() is useless because profiles are already sorted
                    constraint_lines = ConstraintLine.get_lines_and_set_limits_from_sections(
                        section_seq, args.interp_constraint_lines
                    )

                    mesh_constr = MeshConstructor(section_seq=section_seq, lat_step=args.lat_step,
                                                  nb_pts_lat=args.nb_pts_lat, interp_values=args.interp_values)
                    mesh_constr.build_interp(constraint_lines, args.long_step, args.constant_long_disc)
                    mesh_constr.build_mesh(in_floworiented_crs=True)

                    global_mesh_constr.append_mesh_constr(mesh_constr)
                else:
                    logger.warning("Branche ignorée par manque de sections")
            except TatooineException as e:
                logger.error("/!\\ Branche ignorée à cause d'une erreur bloquante :")
                logger.error(e.message)
            logger.info("\n")

    # Handle casiers in floodplain
    nb_casiers = len(modele.get_liste_casiers())
    if args.infile_dem and nb_casiers > 0:
        logger.info("===== TRAITEMENT DES CASIERS =====")

        if not os.path.exists(args.infile_dem):
            raise TatooineException("File not found: %s" % args.infile_dem)
        from gdal import Open
        raster = Open(args.infile_dem)
        dem_interp = interp_raster(raster)

        floodplain_step = args.floodplain_step if not None else args.long_step
        max_elem_area = floodplain_step * floodplain_step / 2.0
        simplify_dist = floodplain_step / 2.0

        for i, casier in enumerate(modele.get_liste_casiers()):
            if casier.is_active:
                if casier.geom is None:
                    raise TatooineException("Geometry of %s could not be found" % casier)
                line = casier.geom.simplify(simplify_dist)
                if not line.is_closed:
                    raise RuntimeError
                coords = resample_2d_line(line.coords, floodplain_step)[1:]  # Ignore last duplicated node

                hard_nodes_xy = np.array(coords, dtype=np.float)
                hard_nodes_idx = np.arange(0, len(hard_nodes_xy), dtype=np.int)
                hard_segments = np.column_stack((hard_nodes_idx, np.roll(hard_nodes_idx, 1)))

                tri = {
                    'vertices': np.array(np.column_stack((hard_nodes_xy[:, 0], hard_nodes_xy[:, 1]))),
                    'segments': hard_segments,
                }
                triangulation = triangle.triangulate(tri, opts='qpa%f' % max_elem_area)

                nodes_xy = np.array(triangulation['vertices'], dtype=np.float)
                bottom = dem_interp(nodes_xy)
                points = unstructured_to_structured(np.column_stack((nodes_xy, bottom)), names=['X', 'Y', 'Z'])

                global_mesh_constr.add_floodplain_mesh(triangulation, points)

    if len(global_mesh_constr.points) == 0:
        raise ExceptionCrue10("Aucun point à traiter, adaptez l'option `--branch_patterns` et/ou `--branch_types_filter`")

    logger.info(global_mesh_constr.summary())  # General information about the merged mesh

    if args.infile_rcal:
        # Read rcal result file
        results = RunResults(args.infile_rcal)
        logger.info(results.summary())

        # Check result consistency
        missing_sections = modele.get_missing_active_sections(results.emh['Section'])
        if missing_sections:
            raise ExceptionCrue10("Sections actives dans le scénario mais manquantes dans le Run :\n%s"
                                  % missing_sections)

        # Subset results to get requested variables at active sections
        varnames_1d = results.variables['Section']
        logger.info("Variables 1D disponibles aux sections: %s" % varnames_1d)
        try:
            pos_z = varnames_1d.index('Z')
        except ValueError:
            raise TatooineException("La variable Z doit être présente dans les résultats aux sections")
        if global_mesh_constr.has_floodplain:
            try:
                pos_z_fp = results.variables['Casier'].index('Z')
            except ValueError:
                raise TatooineException("La variable Z doit être présente dans les résultats aux casiers")
        else:
            pos_z_fp = None

        pos_variables = [results.variables['Section'].index(var) for var in varnames_1d]
        pos_sections_list = [results.emh['Section'].index(profil.id) for profil in global_mesh_constr.section_seq]
        if global_mesh_constr.has_floodplain:
            pos_casiers_list = [results.emh['Casier'].index(casier.id)
                                for casier in modele.get_liste_casiers() if casier.is_active]
        else:
            pos_casiers_list = []

        if 'H' in varnames_1d:
            raise TatooineException("La variable H (charge) de Crue10 entre en conflit avec celle de TatooineMesher. "
                                    "Veuillez supprimer cette variable et relancer le traitement")
        additional_variables_id = ['H']
        if 'Vact' in varnames_1d:
            additional_variables_id.append('M')

        values_geom = global_mesh_constr.interp_values_from_geom()
        z_bottom = values_geom[0, :]
        with Serafin.Write(args.outfile_mesh, args.lang, overwrite=True) as resout:
            title = '%s (written by TatooineMesher)' % os.path.basename(args.outfile_mesh)
            output_header = Serafin.SerafinHeader(title=title, lang=args.lang)
            output_header.from_triangulation(global_mesh_constr.triangle['vertices'],
                                             global_mesh_constr.triangle['triangles'] + 1)
            for var_name in VARIABLES_FROM_GEOMETRY:
                if var_name in ['B', 'W']:
                    output_header.add_variable_from_ID(var_name)
                else:
                    output_header.add_variable_str(var_name, var_name, '')
            for var_id in additional_variables_id:
                output_header.add_variable_from_ID(var_id)
            for var_name in varnames_1d:
                output_header.add_variable_str(var_name, var_name, '')
            resout.write_header(output_header)

            if args.calc_unsteady is None:
                for i, calc_name in enumerate(results.calc_steady_dict.keys()):
                    logger.info("~> Calcul permanent %s" % calc_name)
                    # Read a single *steady* calculation
                    res_steady = results.get_res_steady(calc_name)
                    variables_at_profiles = res_steady['Section'][pos_sections_list, :][:, pos_variables]
                    if global_mesh_constr.has_floodplain:
                        z_at_casiers = res_steady['Casier'][pos_casiers_list, pos_z_fp]
                    else:
                        z_at_casiers = None

                    # Interpolate between sections and set in casiers
                    values_res = global_mesh_constr.interp_values_from_res(variables_at_profiles, z_at_casiers, pos_z)

                    # Compute water depth: H = Z - Zf and clip below 0m (avoid negative values)
                    depth = np.clip(values_res[pos_z, :] - z_bottom, a_min=0.0, a_max=None)

                    # Merge and write values
                    if 'Vact' in varnames_1d:
                        # Compute velocity magnitude from Vact and apply mask "is active bed"
                        velocity = values_res[varnames_1d.index('Vact'), :] * values_geom[1, :]
                        values = np.vstack((values_geom, depth, velocity, values_res))
                    else:
                        values = np.vstack((values_geom, depth, values_res))

                    resout.write_entire_frame(output_header, 3600.0 * i, values)

            else:
                calc_unsteady = results.get_calc_unsteady(args.calc_unsteady)
                logger.info("Calcul transitoire %s" % args.calc_unsteady)
                res_unsteady = results.get_res_unsteady(args.calc_unsteady)

                for i, (time, _) in enumerate(calc_unsteady.frame_list):
                    logger.info("~> %fs" % time)
                    res_at_sections = res_unsteady['Section'][i, :, :]
                    variables_at_profiles = res_at_sections[pos_sections_list, :][:, pos_variables]
                    if global_mesh_constr.has_floodplain:
                        z_at_casiers = res_unsteady['Casier'][i, pos_casiers_list, pos_z_fp]
                    else:
                        z_at_casiers = None

                    # Interpolate between sections
                    values_res = global_mesh_constr.interp_values_from_res(variables_at_profiles, z_at_casiers, pos_z)

                    # Compute water depth: H = Z - Zf and clip below 0m (avoid negative values)
                    depth = np.clip(values_res[pos_z, :] - z_bottom, a_min=0.0, a_max=None)

                    # Merge and write values
                    if 'Vact' in varnames_1d:
                        # Compute velocity magnitude from Vact and apply mask "is active bed"
                        velocity = values_res[varnames_1d.index('Vact'), :] * values_geom[1, :]
                        values = np.vstack((values_geom, depth, velocity, values_res))
                    else:
                        values = np.vstack((values_geom, depth, values_res))

                    resout.write_entire_frame(output_header, time, values)

    else:
        # Write a single frame with only variables from geometry
        global_mesh_constr.export_mesh(args.outfile_mesh, lang=args.lang)

    t2 = perf_counter()
    logger.info("=> Execution time: {}s".format(t2 - t1))


parser = MyArgParse(description=__doc__)
parser.add_common_args(project_straight_line=False, constant_long_disc=True)
# Inputs
parser.infile_args.title = "~> Crue10 input model and run (and the optional DEM)"
parser.infile_args.add_argument("infile_etu", help="Crue10 study file (*.etu.xml)")
parser.infile_args.add_argument("model_name", help="model name")
parser.infile_args.add_argument("--infile_rcal", help="Crue10 results file (*.rcal.xml)")
parser.infile_args.add_argument("--calc_unsteady", help="name of the unsteady file "
                                                        "(otherwise considers all steady calculations)")
parser.infile_args.add_argument("--infile_dem", help='Raster file (geoTIFF format) containing bottom elevation for '
                                                     'the "casiers" in the floodplain (*.tif)')
# Parameters to select branches
parser_branches = parser.add_argument_group("Parameters to filter branches")
parser_branches.add_argument("--branch_types_filter", type=int, nargs='+', default=Branche.TYPES_IN_MINOR_BED,
                             help="types of branches to consider")
parser_branches.add_argument("--branch_patterns", nargs='+', default=None,
                             help="list of patterns to filter branches which name does not contain any pattern")
# Mesh parameters
parser.mesher_args.add_argument("--floodplain_step", type=float, default=None,
                                help="floodplain space step (in m)")
# Outputs
parser.add_out_mesh_file()


if __name__ == '__main__':
    args = parser.parse_args()
    try:
        mesh_crue10_run(args)
    except FileNotFoundError as e:
        logger.critical(e)
        sys.exit(1)
    except ExceptionCrue10 as e:
        logger.critical(e)
        sys.exit(2)
