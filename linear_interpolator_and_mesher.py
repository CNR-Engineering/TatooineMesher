"""
Interpolation linéaire entre profils en travers
Génération d'un maillage (optionnel)

Remarque :
* intégration possible de seuils (correction locale de la bathymétrie)
"""
import time

from core.arg_command_line import MyArgParse
from core.base import LigneContrainte, MeshConstructor, SuiteProfilsTravers
from core.utils import get_axe_hydraulique


def main(args):
    t1 = time.clock()

    print("~> Lecture des fichiers d'entrées")
    axe = get_axe_hydraulique(args.infile_axe)
    profils_travers = SuiteProfilsTravers(args.infile_profils_travers, "Profils en travers", field=args.attr_profils_travers)

    has_epi = args.infile_epis is not None and args.dist_corr_epi is not None
    if has_epi:
        epis = SuiteProfilsTravers(args.infile_epis, "Épi", field=args.attr_epis)

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

    t2 = time.clock()
    print("=> le temps d'execution est de : {}s".format(t2-t1))


if __name__ == '__main__':
    parser = MyArgParse(description=__doc__)
    parser.add_common_args()
    # Inputs
    parser.add_argument("--infile_epis", help="fichier d'entrée i3s/shp avec les épis")
    parser.add_argument("--attr_epis", help="attribute pour identifier les épis")
    parser.add_argument("--dist_corr_epi", type=float,
                        help="distance autour des épis (en général inférieur aux pas trans. et long.)")
    parser.add_argument("--constant_ech_long", help="méthode de calcul du nombre de profils interpolés entre profils "
                        "calculé par profil (constant, ie True) ou par lit (variable, ie False)", action='store_true')
    # Outputs
    parser.add_argument("--outfile_mesh", help="maillage en t3s")
    parser.add_argument("--outfile_semis", help="semis en xyz")
    args = parser.parse_args()

    main(args)
