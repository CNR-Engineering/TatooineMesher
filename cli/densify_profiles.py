#!/usr/bin/env python3
"""
densify_profiles.py

Interpolate intermediate cross-sections
"""
from time import perf_counter

from tatooinemesher.constraint_lines import LigneContrainte
from tatooinemesher.mesh_constructor import MeshConstructor
from tatooinemesher.sections import SuiteProfilsTravers
from tatooinemesher.utils.arg_command_line import MyArgParse
from tatooinemesher.utils import get_axe_hydraulique, logger, set_logger_level, TatooineException


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
        if args.nb_pts_trans is not None and len(lignes_contraintes) != 2:
            raise TatooineException("Argument `--nb_pts_trans` is only compatible with 2 constraint lines!")
    profils_travers.find_and_add_limits(lignes_contraintes, args.dist_max)
    #profils_travers.export_profil_shp('profils_travers_export_profil.shp')

    mesh_constr = MeshConstructor(profils_travers=profils_travers, pas_trans=args.pas_trans,
                                  nb_pts_trans=args.nb_pts_trans)
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
