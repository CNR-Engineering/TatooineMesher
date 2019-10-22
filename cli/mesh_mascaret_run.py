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

from tatooinemesher.constraint_lines import ConstraintLine
from tatooinemesher.mesh_constructor import MeshConstructor
from tatooinemesher.sections import CrossSection, CrossSectionSequence
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
        section_seq = CrossSectionSequence()

        dist_proj_axe = 0.0
        prev_x, prev_y = 0.0, 0.0
        for section_idx, masc_section in enumerate(reach):
            section = CrossSection(masc_section.id,
                                   [(x, y) for x, y in zip(masc_section.x, masc_section.y)],
                                   "Cross-section")

            section.coord.values = np.core.records.fromarrays(
                np.column_stack((masc_section.z,)).T,
                names=VARIABLES_FROM_GEOMETRY
            )
            x, y = masc_section.axis
            if section_idx != 0:
                dist_proj_axe += sqrt((x - prev_x)**2 + (y - prev_y)**2)

            section.dist_proj_axe = dist_proj_axe
            prev_x, prev_y = x, y

            section_seq.add_section(section)

        if len(section_seq) >= 2:
            section_seq.check_intersections()
            # section_seq.sort_by_dist() is useless because cross-sections are already sorted
            constraint_lines = ConstraintLine.get_lines_and_set_limits_from_sections(section_seq,
                                                                                     args.interp_constraint_lines)

            mesh_constr = MeshConstructor(section_seq=section_seq, lat_step=args.lat_step,
                                          nb_pts_lat=args.nb_pts_lat, interp_values=args.interp_values)
            mesh_constr.build_interp(constraint_lines, args.long_step, args.constant_long_disc)
            mesh_constr.build_mesh(in_floworiented_crs=True)

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
                variables_at_sections = masc_res.get_values(idx_time)[reach.id]

                # Interpolate between sections and set in casiers
                values_res = global_mesh_constr.interp_values_from_res(variables_at_sections, None, pos_z)

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
parser.add_common_args(constant_long_disc=True)
# Inputs
parser.infile_args.title = "~> Input Mascaret files arguments"
parser.infile_args.add_argument("infile_geo", help="Mascaret geometry file (*.georef, *.georefC)")
parser.infile_args.add_argument("--infile_res", help="Mascaret results file (*.opt, *.rub)")

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
