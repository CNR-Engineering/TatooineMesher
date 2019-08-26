#!/usr/bin/env python3
"""
Densifier des profils en travers
"""
from time import perf_counter

from core.arg_command_line import MyArgParse
from core.base import LigneContrainte, MeshConstructor, SuiteProfilsTravers
from core.utils import get_axe_hydraulique, logger, set_logger_level, TatooineException


def densify_profiles(args):
    set_logger_level(args.verbose)
    t1 = perf_counter()

    logger.info("~> Lecture des fichiers d'entrées")
    axe = get_axe_hydraulique(args.infile_axe)
    profils_travers = SuiteProfilsTravers.from_file(args.infile_profils_travers, "Profils en travers",
                                                    field_id=args.attr_profils_travers,
                                                    project_straight_line=args.project_straight_line)

    profils_travers.compute_dist_proj_axe(axe, args.dist_max)
    profils_travers.check_intersections()
    profils_travers.sort_by_dist()

    if args.infile_lignes_contraintes is None:
        lignes_contraintes = LigneContrainte.get_lines_and_set_limits_from_profils(profils_travers)
    else:
        lignes_contraintes = LigneContrainte.get_lines_from_file(args.infile_lignes_contraintes)
    profils_travers.find_and_add_limits(lignes_contraintes, args.dist_max)
    #profils_travers.export_profil_shp('profils_travers_export_profil.shp')

    mesh_constr = MeshConstructor(profils_travers, args.pas_trans, args.nb_pts_trans)
    mesh_constr.build_interp(lignes_contraintes, args.pas_long, True)
    mesh_constr.export_profiles(args.outfile_profiles)

    t2 = perf_counter()
    logger.info("=> le temps d'execution est de : {}s".format(t2-t1))


parser = MyArgParse(description=__doc__)
parser.add_common_args()
# Outputs
parser.add_argument("outfile_profiles", help="fichier de sortie contenant les profils en travers (i3s, georefC, shp)")


if __name__ == '__main__':
    args = parser.parse_args()
    try:
        densify_profiles(args)
    except TatooineException as e:
        logger.critical(e.message)
