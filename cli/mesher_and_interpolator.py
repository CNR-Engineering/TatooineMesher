#!/usr/bin/env python3
"""
mesher_and_interpolator.py

Channel mesher and/or interpolator from cross-sections
"""
#TODO: integrate "seuils?" (local correction of the bathymetry)
from time import perf_counter

from tatooinemesher.constraint_lines import LigneContrainte
from tatooinemesher.mesh_constructor import MeshConstructor
from tatooinemesher.sections import ProfilTravers, SuiteProfilsTravers
from tatooinemesher.utils.arg_command_line import MyArgParse
from tatooinemesher.utils import get_axe_hydraulique, logger, set_logger_level, TatooineException


def mesher_and_interpolator(args):
    set_logger_level(args.verbose)
    t1 = perf_counter()

    logger.info("~> Reading input files")
    axe = get_axe_hydraulique(args.infile_axe)
    profils_travers = SuiteProfilsTravers.from_file(args.infile_profils_travers, "Cross-sections",
                                                    field_id=args.attr_profils_travers,
                                                    project_straight_line=args.project_straight_line)

    if args.infile_epis is not None and args.dist_corr_epi is not None:
        epis = SuiteProfilsTravers.from_file(args.infile_epis, "Groynes", field_id=args.attr_epis,
                                             project_straight_line=args.project_straight_line)
    else:
        epis = None

    profils_travers.compute_dist_proj_axe(axe, args.dist_max)
    profils_travers.check_intersections()
    profils_travers.sort_by_dist()
    # profils_travers.export_profil_shp('profiles_projected.shp')  # DEBUG

    if args.infile_lignes_contraintes is None:
        lignes_contraintes = LigneContrainte.get_lines_and_set_limits_from_profils(profils_travers,
                                                                                   args.interp_lignes_contraintes)
    else:
        lignes_contraintes = LigneContrainte.get_lines_from_file(args.infile_lignes_contraintes,
                                                                 args.interp_lignes_contraintes)
        if args.nb_pts_trans is not None and len(lignes_contraintes) != 2:
            raise TatooineException("Argument `--nb_pts_trans` is only compatible with 2 constraint lines!")
        if args.interp_values.startswith('BI') and len(lignes_contraintes) != 2:
            raise TatooineException("A 2D interpolation is only compatible with 2 constraint lines!")
        profils_travers.find_and_add_limits(lignes_contraintes, args.dist_max)

    mesh_constr = MeshConstructor(profils_travers=profils_travers, pas_trans=args.pas_trans,
                                  nb_pts_trans=args.nb_pts_trans, interp_values=args.interp_values)
    mesh_constr.build_interp(lignes_contraintes, args.pas_long, args.constant_ech_long)
    # mesh_constr.export_segments('check_segments.shp')  # DEBUG

    if epis is not None:
        mesh_constr.corr_bathy_on_epis(epis, args.dist_corr_epi)

    if args.outfile_semis is not None:
        mesh_constr.export_points(args.outfile_semis)

    if args.outfile_mesh is not None:
        mesh_constr.build_mesh()
        mesh_constr.export_mesh(args.outfile_mesh)

    t2 = perf_counter()
    logger.info("=> Execution time: {}s".format(t2-t1))


parser = MyArgParse(description=__doc__)
parser.add_common_args()
# Inputs
parser_epis = parser.add_argument_group('Parameters to define lateral groynes (optional)')
parser_epis.add_argument("--infile_epis", help="input file for groynes i3s/shp")
parser_epis.add_argument("--attr_epis", help="attribute to identify groynes")
parser_epis.add_argument("--dist_corr_epi", type=float,
                         help="distance around groynes to modify nodes close to them "
                              "(should be less than lateral and longitudinal space step)")
parser.add_argument("--constant_ech_long", help="méthode de calcul du nombre de profils interpolés entre profils : "
                    "par profil (constant, ie True) ou par lit (variable, ie False)", action='store_true')
parser.add_argument("--interp_lignes_contraintes",
                    help="méthode d'interpolation des coordonnées des lignes de contraintes",
                    default='LINEAR', choices=('LINEAR', 'FINITE_DIFF', 'CARDINAL'))
parser.add_argument("--interp_values",
                    help="méthode d'interpolation (transversale ou 2D) pour les valeurs",
                    default='LINEAR', choices=('LINEAR', 'B-SPLINE', 'AKIMA', 'PCHIP', 'CUBIC_SPLINE',
                                               'BILINEAR', 'BICUBIC', 'BIVARIATE_SPLINE'))
# Outputs
parser_outfiles = parser.add_argument_group('Output files')
parser_outfiles.add_argument("--outfile_mesh", help="mesh file. Supported formats: slf (Telemac), t3s (BlueKenue),"
                                                    " and xml (LandXML for GIS)")
parser_outfiles.add_argument("--outfile_semis", help="semis de points au format xyz"
                                                     " (correspond to mesh nodes)")


if __name__ == '__main__':
    args = parser.parse_args()
    try:
        mesher_and_interpolator(args)
    except TatooineException as e:
        logger.critical(e.message)
