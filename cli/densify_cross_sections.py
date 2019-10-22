#!/usr/bin/env python3
"""
densify_cross_sections.py

Interpolate initial and intermediate cross-sections
Multiple variables are supported if input cross-sections file is a shapefile with POINTZ type.
"""
from time import perf_counter

from tatooinemesher.constraint_line import ConstraintLine
from tatooinemesher.mesh_constructor import MeshConstructor
from tatooinemesher.section import CrossSectionSequence
from tatooinemesher.utils.arg_command_line import MyArgParse
from tatooinemesher.utils import get_hydraulic_axis, logger, set_logger_level, TatooineException


def densify_cross_sections(args):
    set_logger_level(args.verbose)
    t1 = perf_counter()

    logger.info("~> Reading input files")
    axe = get_hydraulic_axis(args.infile_axis)
    section_seq = CrossSectionSequence.from_file(args.infile_cross_sections, "Cross-section",
                                                 field_id=args.attr_cross_sections,
                                                 project_straight_line=args.project_straight_line)

    section_seq.compute_dist_proj_axe(axe, args.dist_max)
    section_seq.check_intersections()
    section_seq.sort_by_dist()

    if args.infile_constraint_lines is None:
        constraint_lines = ConstraintLine.get_lines_and_set_limits_from_sections(section_seq,
                                                                                 args.interp_constraint_lines)
    else:
        constraint_lines = ConstraintLine.get_lines_from_file(args.infile_constraint_lines,
                                                              args.interp_constraint_lines)
        if args.nb_pts_lat is not None and len(constraint_lines) != 2:
            raise TatooineException("Argument `--nb_pts_lat` is only compatible with 2 constraint lines!")
        if args.interp_values.startswith('BI') and len(constraint_lines) != 2:
            raise TatooineException("A 2D interpolation is only compatible with 2 constraint lines!")
        section_seq.find_and_add_limits(constraint_lines, args.dist_max)

    #section_seq.export_sections_shp('export_cross-sections.shp')

    mesh_constr = MeshConstructor(section_seq=section_seq, lat_step=args.lat_step,
                                  nb_pts_lat=args.nb_pts_lat, interp_values=args.interp_values)
    mesh_constr.build_interp(constraint_lines, args.long_step, True)
    mesh_constr.export_sections(args.outfile_sections)

    t2 = perf_counter()
    logger.info("=> Execution time: {}s".format(t2 - t1))


parser = MyArgParse(description=__doc__)
parser.add_common_args(project_straight_line=True)
# Inputs
parser.infile_args.add_argument("infile_axis", help="hydraulic axis file (*.shp, *.i2s)")
parser.infile_args.add_argument("infile_cross_sections",
                                help="cross-sections file (*.shp, *.i3s)")
parser.infile_args.add_argument("--infile_constraint_lines",
                                help="constraint lines file (*.shp, *.i2s)")
parser.infile_args.add_argument("--attr_cross_sections", help="attribute to identify cross-sections")
# Outputs
parser.outfile_args.add_argument("outfile_sections",
                                 help="output file containing interpolated cross-sections (*.i3s, *.georefC, *.shp)")


if __name__ == '__main__':
    args = parser.parse_args()
    try:
        densify_cross_sections(args)
    except TatooineException as e:
        logger.critical(e.message)
