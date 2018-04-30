"""

* intégration de seuils (correction locale de la bathymétrie)

Remarque : les abscisses recalculées par l'outil sont identiques à celle données par le module shapely (pour les LineString).

Code source repris de :
J:\DI-Affaires-2014\I.00846.001 - Endiguement chambéry\4- déroulement d'affaire\4-3 mission\5. Hydraulique\01-Domaine&topo\lit_mineur\_scripts avec l'exemple de la Leysse amont
"""
from jinja2 import Environment, FileSystemLoader
from math import ceil
import numpy as np
import os.path
from shapely.geometry import Point
import sys
import time
import triangle

from core.arg_command_line import MyArgParse
from core.base import SuiteProfilsTravers, LigneContrainte, Lit, MeshConstructor
from core.utils import get_axe_hydraulique


# PARAMETRES
DIGITS = 4
SEP = ';'


def main(args):
    t1 = time.clock()

    print("~> Lecture des fichiers d'entrées")
    axe = get_axe_hydraulique(args.infile_axe)
    profils_travers_ori = SuiteProfilsTravers(args.infile_profils_travers, "Profils en travers", field=args.attr_profils_travers)
    lignes_contraintes = LigneContrainte.get_lines_from_file(args.infile_lignes_contraintes)

    has_epi = args.infile_epis is not None and args.dist_corr_epi is not None
    if has_epi:
        epis = SuiteProfilsTravers(args.infile_epis, "Épi", field=args.attr_epis)

    print("~> Calcul des abscisses sur l'axe hydraulique (pour ordonner les profils/épis)")
    profils_travers_ori.calculer_dist_proj_axe(axe)
    profils_dist_proj_axe = np.array([profil.dist_proj_axe for profil in profils_travers_ori])

    print("~> Recherche des limites de lits")
    profils_travers_ori.find_and_add_limits(lignes_contraintes, args.dist_max)

    print("~> Vérifications non intersections des profils et des épis")
    profils_travers = profils_travers_ori
    profils_travers.check_intersections()

    print("~> Classement des profils (et épis) le long de l'axe hydraulique")
    profils_travers.calculer_dist_proj_axe(axe)
    profils_travers = sorted(profils_travers, key=lambda x: x.dist_proj_axe)

    mesh_constr = MeshConstructor()
    mesh_constr.interp(profils_travers, lignes_contraintes, args.pas_trans, args.pas_long, args.constant_ech_long)

    if has_epi:
        mesh_constr.corr_epis(epis, args.dist_corr_epi)

    print("~> Exports en xyz puis en CSV des données")

    with open(args.outfile_semis, 'wb') as fileout:
        np.savetxt(fileout, mesh_constr.points, fmt='%.4f')

    first_profil = True
    for profil in profils_travers:
        profil.export_profil_csv(args.outfile_profils, first_profil, SEP, DIGITS)
        profil.export_limites_csv(args.outfile_limites, first_profil, SEP, DIGITS)
        if first_profil:
            first_profil = False

    print("~> Calcul du maillage")

    tri = mesh_constr.export_as_dict()
    t = triangle.triangulate(tri, opts='p')
    nnode = len(t['vertices'])
    nelem = len(t['triangles'])
    print("Génération d'un maillage avec {} noeuds et {} éléments".format(nnode, nelem))

    print("~> Écriture du maillage")

    out_extension = args.outfile_mesh[-4:]
    if out_extension == ".t3s":
        with open(args.outfile_mesh, 'w', newline='') as fileout:
            # Écriture en-tête
            date = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
            fileout.write("""#########################################################################
    :FileType t3s  ASCII  EnSim 1.0
    # Canadian Hydraulics Centre/National Research Council (c) 1998-2012
    # DataType                 2D T3 Scalar Mesh
    #
    :Application              BlueKenue
    :Version                  3.3.4
    :WrittenBy                scripts_LDN
    :CreationDate             {}
    #
    #------------------------------------------------------------------------
    #
    :NodeCount {}
    :ElementCount {}
    :ElementType  T3
    #
    :EndHeader
    """.format(date, nnode, nelem))

        with open(args.outfile_mesh, mode='ab') as fileout:
            # Tableau des coordonnées (x, y, z)
            np.savetxt(fileout, np.column_stack((t['vertices'], mesh_constr.points['Z'])), delimiter=' ', fmt='%.{}f'.format(DIGITS))

            # Tableau des éléments (connectivité)
            np.savetxt(fileout, t['triangles']+1, delimiter=' ', fmt='%i')

    elif out_extension == ".xml":
        env = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')))
        template = env.get_template("LandXML_template.xml")
        template_render = template.render(
            nodes=np.round(np.column_stack((t['vertices'], mesh_constr.points['Z'])), DIGITS),
            ikle=t['triangles']+1
        )

        # Écriture du fichier XML
        with open(args.outfile_mesh, 'w') as fileout:
            fileout.write(template_render)

    t2 = time.clock()
    print("=> le temps d'execution est de : {}s".format(t2-t1))


if __name__ == "__main__":
    parser = MyArgParse(description=__doc__)
    parser.add_common_args()
    # Inputs
    parser.add_argument("--infile_epis", help="fichier d'entrée i3s/shp avec les épis")
    parser.add_argument("--attr_epis", help="attribute pour identifier les épis")
    parser.add_argument("--dist_corr_epi", type=float,
                        help="distance autour des épis (en général inférieur aux pas trans. et long.)")
    # Outputs
    parser.add_argument("outfile_mesh", help="maillage en t3s")
    parser.add_argument("outfile_semis", help="semis en xyz")
    parser.add_argument("outfile_limites", help="limites en csv")
    parser.add_argument("outfile_profils", help="profils en csv")
    args = parser.parse_args()

    main(args)
