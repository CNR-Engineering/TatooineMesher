#!/usr/bin/env python3
"""
Interpolation linéaire entre profils en travers
Génération d'un maillage (optionnel)

Remarque :
* intégration possible de seuils (correction locale de la bathymétrie)
"""
import time

from core.arg_command_line import MyArgParse
from core.base import LigneContrainte, MeshConstructor, SuiteProfilsTravers
from core.utils import get_axe_hydraulique, logger, set_logger_level, TatooineException


def linear_interpolator_and_mesher(args):
    set_logger_level(args.verbose)
    t1 = time.process_time()

    logger.info("~> Lecture des fichiers d'entrées")
    axe = get_axe_hydraulique(args.infile_axe)
    profils_travers = SuiteProfilsTravers.from_file(args.infile_profils_travers, "Profils en travers",
                                                    field_id=args.attr_profils_travers)

    has_epi = args.infile_epis is not None and args.dist_corr_epi is not None
    if has_epi:
        epis = SuiteProfilsTravers.from_file(args.infile_epis, "Épi", field_id=args.attr_epis)

    profils_travers.compute_dist_proj_axe(axe)
    profils_travers.check_intersections()
    profils_travers.sort_by_dist()

    if args.infile_lignes_contraintes is None:
        lignes_contraintes = LigneContrainte.get_lines_from_profils(profils_travers)
    else:
        lignes_contraintes = LigneContrainte.get_lines_from_file(args.infile_lignes_contraintes)
    profils_travers.find_and_add_limits(lignes_contraintes, args.dist_max)

    mesh_constr = MeshConstructor(profils_travers, args.pas_trans)
    mesh_constr.build_interp(lignes_contraintes, args.pas_long, args.constant_ech_long)

    if has_epi:
        mesh_constr.corr_bathy_on_epis(epis, args.dist_corr_epi)

    if args.outfile_semis is not None:
        mesh_constr.export_points(args.outfile_semis)

    if args.outfile_mesh is not None:
        mesh_constr.build_mesh()
        mesh_constr.export_mesh(args.outfile_mesh)

    t2 = time.process_time()
    logger.info("=> le temps d'execution est de : {}s".format(t2-t1))


parser = MyArgParse(description=__doc__)
parser.add_common_args()
# Inputs
parser_epis = parser.add_argument_group('Paramètres pour la définition des épis (optionnels)')
parser_epis.add_argument("--infile_epis", help="fichier d'entrée i3s/shp avec les épis")
parser_epis.add_argument("--attr_epis", help="attribut pour identifier les épis")
parser_epis.add_argument("--dist_corr_epi", type=float,
                         help="distance autour des épis (en général inférieur aux pas trans. et long.)")
parser.add_argument("--constant_ech_long", help="méthode de calcul du nombre de profils interpolés entre profils : "
                    "par profil (constant, ie True) ou par lit (variable, ie False)", action='store_true')
# Outputs
parser_outfiles = parser.add_argument_group('Fichiers de sortie')
parser_outfiles.add_argument("--outfile_mesh", help="maillage au format : slf (Telemac), t3s (BlueKenue),"
                                                    " ou xml (LandXML pour ArcGIS)")
parser_outfiles.add_argument("--outfile_semis", help="semis de points au format xyz"
                                                     " (correspond aux noeuds du maillage")


if __name__ == '__main__':
    args = parser.parse_args()
    try:
        linear_interpolator_and_mesher(args)
    except TatooineException as e:
        logger.critical(e.message)
