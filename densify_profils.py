"""
Densifier des profils en travers
"""
import time

from core.arg_command_line import MyArgParse
from core.utils import get_axe_hydraulique
from core.base import LigneContrainte, MeshConstructor, SuiteProfilsTravers


def main(args):
    t1 = time.clock()

    print("~> Lecture des fichiers d'entrées")
    axe = get_axe_hydraulique(args.infile_axe)
    profils_travers = SuiteProfilsTravers(args.infile_profils_travers, "Profils en travers", field_id=args.attr_profils_travers)

    profils_travers.compute_dist_proj_axe(axe)
    profils_travers.check_intersections()
    profils_travers.sort_by_dist()

    if args.infile_lignes_contraintes is None:
        lignes_contraintes = LigneContrainte.get_lines_from_profils(profils_travers)
    else:
        lignes_contraintes = LigneContrainte.get_lines_from_file(args.infile_lignes_contraintes)
    profils_travers.find_and_add_limits(lignes_contraintes, args.dist_max)
    #profils_travers.export_profil_shp('profils_travers_export_profil.shp')

    mesh_constr = MeshConstructor(profils_travers, args.pas_trans)
    mesh_constr.build_interp(lignes_contraintes, args.pas_long, True)
    mesh_constr.export_profiles(args.out_profiles)

    t2 = time.clock()
    print("=> le temps d'execution est de : {}s".format(t2-t1))


if __name__ == '__main__':
    parser = MyArgParse(description=__doc__)
    parser.add_common_args()
    args = parser.parse_args()
    # Outputs
    parser.add_argument("outfile_profiles", help="fichier de sortie contenant les profils en travers (i3s)")

    main(args)
