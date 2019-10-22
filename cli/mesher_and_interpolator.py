#!/usr/bin/env python3
"""
mesher_and_interpolator.py

Channel mesher and/or interpolator from cross-sections
"""
#TODO: integrate "seuils?" (local correction of the bathymetry)
from time import perf_counter

from tatooinemesher.constraint_lines import ConstraintLine
from tatooinemesher.mesh_constructor import MeshConstructor
from tatooinemesher.sections import CrossSection, CrossSectionSequence
from tatooinemesher.utils.arg_command_line import MyArgParse
from tatooinemesher.utils import get_hydraulic_axis, logger, set_logger_level, TatooineException


def mesher_and_interpolator(args):
    set_logger_level(args.verbose)
    t1 = perf_counter()

    logger.info("~> Reading input files")
    axe = get_hydraulic_axis(args.infile_axis)
    section_seq = CrossSectionSequence.from_file(args.infile_cross_sections, "Cross-section",
                                                 field_id=args.attr_cross_sections,
                                                 project_straight_line=args.project_straight_line)

    # if args.infile_epis is not None and args.dist_corr_epi is not None:
    #     epis = CrossSectionSequence.from_file(args.infile_epis, "Groynes", field_id=args.attr_epis,
    #                                           project_straight_line=args.project_straight_line)
    # else:
    #     epis = None

    section_seq.compute_dist_proj_axe(axe, args.dist_max)
    section_seq.check_intersections()
    section_seq.sort_by_dist()
    # section_seq.export_sections_shp('profiles_projected.shp')  # DEBUG

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

    mesh_constr = MeshConstructor(section_seq=section_seq, lat_step=args.lat_step,
                                  nb_pts_lat=args.nb_pts_lat, interp_values=args.interp_values)
    mesh_constr.build_interp(constraint_lines, args.long_step, args.constant_long_disc)
    # mesh_constr.export_segments('check_segments.shp')  # DEBUG

    # if epis is not None:
    #     mesh_constr.corr_bathy_on_epis(epis, args.dist_corr_epi)

    if args.outfile_nodes is not None:
        mesh_constr.export_points(args.outfile_nodes)

    if args.outfile_mesh is not None:
        mesh_constr.build_mesh()
        mesh_constr.export_mesh(args.outfile_mesh)

    t2 = perf_counter()
    logger.info("=> Execution time: {}s".format(t2-t1))


parser = MyArgParse(description=__doc__)
parser.add_common_args(project_straight_line=True, constant_long_disc=True)
# Inputs
parser.infile_args.add_argument("infile_axis", help="hydraulic axis file (*.shp, *.i2s)")
parser.infile_args.add_argument("infile_cross_sections",
                                help="cross-sections file (*.shp, *.i3s)")
parser.infile_args.add_argument("--infile_constraint_lines",
                                help="constraint lines file (*.shp, *.i2s)")
parser.infile_args.add_argument("--attr_cross_sections", help="attribute to identify cross-sections")
#TODO: add groynes
# parser_epis = parser.add_argument_group('Parameters to define lateral groynes (optional)')
# parser_epis.add_argument("--infile_epis", help="input file for groynes (*.shp, *.i3s)")
# parser_epis.add_argument("--attr_epis", help="attribute to identify groynes")
# parser_epis.add_argument("--dist_corr_epi", type=float,
#                          help="distance around groynes to modify nodes close to them "
#                               "(should be less than lateral and longitudinal space step)")
# Outputs
parser_outfiles = parser.add_argument_group('Output files')
parser_outfiles.add_argument("--outfile_mesh",
                             help="output mesh file. Supported formats: slf (Telemac), t3s (BlueKenue), "
                                  "and xml (LandXML for GIS)")
parser_outfiles.add_argument("--outfile_nodes", help="output points set file with mesh nodes (*.shp, *.xyz)")


if __name__ == '__main__':
    args = parser.parse_args()
    try:
        mesher_and_interpolator(args)
    except TatooineException as e:
        logger.critical(e.message)
