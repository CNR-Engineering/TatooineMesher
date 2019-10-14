"""
mesh_mascaret_run.py

Mesh Mascaret geometry and results if specified
"""
from math import sqrt
import numpy as np
import os.path
from pyteltools.slf import Serafin
import sys
from time import perf_counter

from crue10.utils import CrueError
from mascaret.mascaret_file import MascaretFile
from mascaret.mascaretgeo_file import MascaretGeoFile

from tatooinemesher.constraint_lines import LigneContrainte
from tatooinemesher.mesh_constructor import MeshConstructor
from tatooinemesher.sections import ProfilTravers, SuiteProfilsTravers
from tatooinemesher.utils.arg_command_line import MyArgParse
from tatooinemesher.utils import logger, set_logger_level, TatooineException


VARIABLES_FROM_GEOMETRY = ['B']


def mesh_mascaret_run(args):
    set_logger_level(args.verbose)
    t1 = perf_counter()

    masc_geo = MascaretGeoFile(args.infile_geo)
    logger.info("Read %s " % masc_geo)
    # masc_geo.export_shp_lines(args.infile_geo.replace('.georef', '.shp'))
    if not masc_geo.has_ref:
        raise TatooineException("The file `%s` does not contain any georeferenced data" % masc_geo.file_name)

    global_mesh_constr = MeshConstructor()

    for reach_id, reach in masc_geo.reaches.items():
        logger.info(reach)
        profils_travers = SuiteProfilsTravers()

        dist_proj_axe = 0.0
        prev_x, prev_y = 0.0, 0.0
        for idx_section, section in enumerate(reach):
            profile = ProfilTravers(section.id, [(x, y) for x, y in zip(section.x, section.y)], 'Section')

            profile.coord.values = np.core.records.fromarrays(
                np.column_stack((section.z,)).T,
                names=VARIABLES_FROM_GEOMETRY
            )
            x, y = section.axis
            if idx_section != 0:
                dist_proj_axe += sqrt((x - prev_x)**2 + (y - prev_y)**2)

            profile.dist_proj_axe = dist_proj_axe
            prev_x, prev_y = x, y

            profils_travers.add_profile(profile)

        if len(profils_travers) >= 2:
            profils_travers.check_intersections()
            # profils_travers.sort_by_dist() is useless because profiles are already sorted
            lignes_contraintes = LigneContrainte.get_lines_and_set_limits_from_profils(profils_travers)

            mesh_constr = MeshConstructor(profils_travers=profils_travers,
                                          pas_trans=args.pas_trans, nb_pts_trans=args.nb_pts_trans)
            mesh_constr.build_interp(lignes_contraintes, args.pas_long, args.constant_ech_long)
            mesh_constr.build_mesh(True)

            global_mesh_constr.append_mesh_constr(mesh_constr)
        else:
            logger.error("/!\\ Reach %s ignored because it does not contain at least 2 sections" % reach_id)

    if len(global_mesh_constr.points) == 0:
        raise CrueError("No node in the generated mesh!")

    logger.info(global_mesh_constr.summary())  # General information about the merged mesh

    if args.infile_res:
        masc_res = MascaretFile(args.infile_res)
        masc_res.get_reaches()
        nb_section_in_geom = masc_geo.nsections
        if masc_res.nsections != nb_section_in_geom:
            raise TatooineException("The number of sections is different between geometry (%i) and results file (%i)"
                                    % (nb_section_in_geom, masc_res.nsections))

        varnames_1d = masc_res.varnames_dict['abbr']
        logger.info("Variables 1D available at sections: %s" % varnames_1d)
        try:
            pos_z = varnames_1d.index('Z')
        except ValueError:
            raise TatooineException("The variable Z must be present in the results file")

        additional_variables_id = ['H']

        values_geom = global_mesh_constr.interp_values_from_geom()
        z_bottom = values_geom[0, :]
        with Serafin.Write(args.outfile_mesh, args.lang, overwrite=True) as resout:
            title = '%s (written by TatooineMesher)' % os.path.basename(args.outfile_mesh)
            output_header = Serafin.SerafinHeader(title=title, lang=args.lang)
            output_header.from_triangulation(global_mesh_constr.triangle['vertices'],
                                             global_mesh_constr.triangle['triangles'] + 1)
            for var_name in VARIABLES_FROM_GEOMETRY:
                if var_name == 'B':
                    output_header.add_variable_from_ID(var_name)
                else:
                    output_header.add_variable_str(var_name, var_name, '')
            for var_id in additional_variables_id:
                output_header.add_variable_from_ID(var_id)
            for var_name in varnames_1d:
                output_header.add_variable_str(var_name, var_name, '')
            resout.write_header(output_header)

            for idx_time, time in enumerate(masc_res.times):
                variables_at_profiles = masc_res.get_values(idx_time)[reach.id]

                # Interpolate between sections and set in casiers
                values_res = global_mesh_constr.interp_values_from_res(variables_at_profiles, None, pos_z)

                # Compute water depth: H = Z - Zf and clip below 0m (avoid negative values)
                depth = np.clip(values_res[pos_z, :] - z_bottom, a_min=0.0, a_max=None)

                values = np.vstack((values_geom, depth, values_res))
                resout.write_entire_frame(output_header, time, values)

    else:
        # Write a single frame with only variables from geometry
        global_mesh_constr.export_mesh(args.outfile_mesh, lang=args.lang)

    t2 = perf_counter()
    logger.info("=> Execution time: {}s".format(t2 - t1))


parser = MyArgParse(description=__doc__)
# Inputs
parser_infiles = parser.add_argument_group("Mascaret input files")
parser_infiles.add_argument("infile_geo", help="Mascaret geometry file (*.georef, *.georefC)")
parser_infiles.add_argument("--infile_res", help="Mascaret results file (*.opt, *.rub)")

# Mesh parameters
parser_mesh = parser.add_argument_group("Parameters to generate the 2D mesh")
parser_mesh.add_argument("--dist_max", type=float, help="maximum search distance to rescue intersections "
                                                        "for limits (in m)", default=0.01)
parser.add_argument("--pas_long", type=float, help="longitudinal space step (in m)")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--pas_trans", type=float, help="lateral space step (in m)")
group.add_argument("--nb_pts_trans", type=int, help="number of nodes laterally")
parser_mesh.add_argument("--constant_ech_long",
                         help="method to compute the number of intermediate profiles: "
                              "per profile (constant, ie True) or per bed (variable, ie False)",
                         action='store_true')

# Outputs
parser_outfiles = parser.add_argument_group('Output file')
parser_outfiles.add_argument("outfile_mesh", help="Telemac results file")
parser_outfiles.add_argument("--lang", help="Language for standard variables in output file", default='fr',
                             choices=['en', 'fr'])

if __name__ == '__main__':
    args = parser.parse_args()
    try:
        mesh_mascaret_run(args)
    except FileNotFoundError as e:
        logger.critical(e)
        sys.exit(1)
    except CrueError as e:
        logger.critical(e)
        sys.exit(2)
