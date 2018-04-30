# -*- coding: utf-8 -*-
"""

* intégration de seuils (correction locale de la bathymétrie)

Remarque : les abscisses recalculées par l'outil sont identiques à celle données par le module shapely (pour les LineString).

Code source repris de :
J:\DI-Affaires-2014\I.00846.001 - Endiguement chambéry\4- déroulement d'affaire\4-3 mission\5. Hydraulique\01-Domaine&topo\lit_mineur\_scripts avec l'exemple de la Leysse amont
"""

# Compatibility with Python2
from __future__ import print_function

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
sep = ';'

parser = MyArgParse(description=__doc__)
# Input
parser.add_argument("infile_axe", help="fichier d'entrée i2s de l'axe hydraulique")
parser.add_argument("infile_profils_travers", help="fichier d'entrée i3s de profils en travers")
parser.add_argument("infile_lignes_contraintes", help="fichier d'entrée i2s de lignes de contrainte")
parser.add_argument("--infile_epis", help="fichier d'entrée i3s avec les épis")
parser.add_argument("--dist_corr_epi", type=float, help="distance autour des épis (en général inférieur aux pas trans. et long.)")
parser.add_argument("--dist_max", type=float, help="distance de recherche maxi des 'intersections fictifs' pour les limites de lits", default=0.01)
parser.add_argument("--constant_ech_long", help="méthode de calcul du nombre de profils interpolés entre profils calculé par profil (constant, ie True) ou par lit (variable, ie False)", action='store_true')
parser.add_argument("--pas_long", type=float, help="pas d'interpolation longitudinal en m.", default=5)
parser.add_argument("--pas_trans", type=float, help="pas d'interpolation transversal en m.", default=3.5)

# args.constant_ech_long :
# - True = nombre de profils interpolés calculé sur les distances minimales entre profils
# - False = nombre des profils interpolés calculé pour chaque lit sur la distance minimale

# Output
parser.add_argument("outfile_mesh", help="maillage en t3s")
parser.add_argument("outfile_semis", help="semis en xyz")
parser.add_argument("outfile_limites", help="limites en csv")
parser.add_argument("outfile_profils", help="profils en csv")

args = parser.parse_args()


# MAIN
t1 = time.clock()

print("~> Lecture des fichiers d'entrées")

# Parcours des fichiers d'entrée
axe = get_axe_hydraulique(args.infile_axe)
profils_travers_ori = SuiteProfilsTravers(args.infile_profils_travers, "Profils en travers", value_as_id=True)
lignes_contraintes = LigneContrainte.get_lines_from_file(args.infile_lignes_contraintes)

if args.infile_epis is not None and args.dist_corr_epi is not None:
    has_epi = True
else:
    has_epi = False

if has_epi:
    epis_ori = SuiteProfilsTravers(args.infile_epis, "Épi", value_as_id=True)

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

print("~> Interpolation sur les profils existants en prenant en compte le passage des lignes de contraintes")
for i in range(len(profils_travers)):
    cur_profil = profils_travers[i]
    print(cur_profil)

    # Recherche des limites communes "amont"
    if i == 0:
        common_limites_id_1 = cur_profil.limites.index
    else:
        common_limites_id_1 = cur_profil.common_limits(profils_travers[i-1])

    # Recherche des limites communes "aval"
    if i == len(profils_travers)-1:
        common_limites_id_2 = cur_profil.limites.index
    else:
        common_limites_id_2 = cur_profil.common_limits(profils_travers[i+1])

    # Union (non ordonnée) des limites amont/aval
    limites_id = np.union1d(common_limites_id_1, common_limites_id_2)
    # Ré-ordonne les limites
    limites_id = cur_profil.limites.index[np.in1d(cur_profil.limites.index, limites_id, assume_unique=True)]

    first_lit = True
    for id1, id2 in zip(limites_id, limites_id[1:]):
        lit = cur_profil.extraire_lit(id1, id2)
        coord_int = lit.interp_along_lit_auto(args.pas_trans)

        if first_lit:
            cur_profil.limites.loc[id1,'id_pt'] = mesh_constr.i_pt+1
        else:
            coord_int = coord_int[1:]

        mesh_constr.add_points(coord_int)

        cur_profil.limites.loc[id2,'id_pt'] = mesh_constr.i_pt

        # Ajoute les nouveaux segments
        new_i_pt = np.arange(cur_profil.limites['id_pt'].loc[id1],
                             cur_profil.limites['id_pt'].loc[id2]+1)
        mesh_constr.add_segments_from_node_list(new_i_pt)

        if first_lit:
            first_lit = False


first_profil = True
### BOUCLE SUR L'ESPACE INTER-PROFIL
print("~> Construction du maillage par zone interprofils puis par lit")
for i, (prev_profil, next_profil) in enumerate(zip(profils_travers, profils_travers[1:])):
    print("> Zone n°{} : entre {} et {}".format(i, prev_profil, next_profil))

    if args.constant_ech_long:
        nb_pts_inter = prev_profil.calcul_nb_pts_inter(next_profil, args.pas_long)
        Xp_adm_list = np.linspace(0.0, 1.0, num=nb_pts_inter+2)[1:-1]

    # Recherche des limites communes entre les deux profils
    common_limites_id = prev_profil.common_limits(next_profil)
    print("Limites de lits communes : {}".format(list(common_limites_id)))

    if len(common_limites_id) < 2:
        sys.exit("ERREUR: aucune interpolation pour l'intervalle {} ({} limites communes)".format(i, len(common_limites_id)))

    else:
        first_lit = True
        ### BOUCLE SUR LES LITS (= MORCEAU(X) DE PROFIL)
        for id1, id2 in zip(common_limites_id, common_limites_id[1:]):
            pt_list_L1 = []
            pt_list_L2 = []

            print("Lit {}-{}".format(id1, id2))

            # Extraction d'une partie des profils
            lit_1 = prev_profil.extraire_lit(id1, id2)
            lit_2 = next_profil.extraire_lit(id1, id2)

            # Abscisses curvilignes le long de la ligne de contrainte
            (Xp_profil1_L1, Xp_profil1_L2) = prev_profil.get_Xt_lignes(id1, id2)
            (Xp_profil2_L1, Xp_profil2_L2) = next_profil.get_Xt_lignes(id1, id2)
            dXp_L1 = Xp_profil2_L1 - Xp_profil1_L1
            dXp_L2 = Xp_profil2_L2 - Xp_profil1_L2

            if dXp_L1 < 0:
                sys.exit("La ligne {} n'est pas orientée dans le même ordre que les profils".format(id1))
            if dXp_L2 < 0:
                sys.exit("La ligne {} n'est pas orientée dans le même ordre que les profils".format(id2))

            if not args.constant_ech_long:
                nb_pts_inter = ceil(min(dXp_L1, dXp_L2)/args.pas_long) - 1
                Xp_adm_list = np.linspace(0.0, 1.0, num=nb_pts_inter+2)[1:-1]

            L1_coord_int = lignes_contraintes[id1].coord_sampling_along_line(Xp_profil1_L1, Xp_profil2_L1, Xp_adm_list)
            L2_coord_int = lignes_contraintes[id2].coord_sampling_along_line(Xp_profil1_L2, Xp_profil2_L2, Xp_adm_list)

            ### BOUCLE SUR LES LIGNES
            for j in range(nb_pts_inter):
                Xp = Xp_adm_list[j]
                P1 = Point(tuple(L1_coord_int[j]))
                P2 = Point(tuple(L2_coord_int[j]))

                lit_int = Lit(lit_1.interp_inter_lineaire(lit_2, Xp, ceil(P1.distance(P2)/args.pas_trans)+1), ['Xt', 'xt'])
                lit_int.move_between_targets(P1, P2)
                coord_int = lit_int.array[['X', 'Y', 'Z']]
                pt_list_L1.append(mesh_constr.i_pt+1)

                if not first_lit:
                    # ignore le 1er point car la ligne de contrainte a déjà été traitée
                    coord_int = coord_int[1:]

                mesh_constr.add_points(coord_int)

                pt_list_L2.append(mesh_constr.i_pt)

            pt_list_L2 = np.array([prev_profil.limites['id_pt'].loc[id2]] + pt_list_L2 + [next_profil.limites['id_pt'].loc[id2]])
            mesh_constr.add_segments_from_node_list(pt_list_L2)

            if first_lit:
                pt_list_L1 = np.array([prev_profil.limites['id_pt'].loc[id1]] + pt_list_L1 + [next_profil.limites['id_pt'].loc[id1]])
                mesh_constr.add_segments_from_node_list(pt_list_L1)
                first_lit = False

    if first_profil:
        first_profil = False

print("~> Correction de la bathymétrie autour des épis")
if has_epi:
    for epi in epis_ori:
        epi_geom = epi.coord.convert_as_linestring()
        for i, coord in enumerate(mesh_constr.points):
            pt_node = Point(tuple(coord))
            if epi_geom.distance(pt_node) < args.dist_corr_epi:
                Xt_proj = epi_geom.project(pt_node)
                pt_proj = epi_geom.interpolate(Xt_proj)
                mesh_constr.points['Z'][i] = pt_proj.z

print("~> Exports en xyz puis en CSV des données")

with open(args.outfile_semis, 'wb') as fileout:
    np.savetxt(fileout, mesh_constr.points, fmt='%.4f')

first_profil = True
for profil in profils_travers:
    profil.export_profil_csv(args.outfile_profils, first_profil, sep, DIGITS)
    profil.export_limites_csv(args.outfile_limites, first_profil, sep, DIGITS)
    if first_profil:
        first_profil = False

print("~> Calcul du maillage")

tri = mesh_constr.export_as_dict()
t = triangle.triangulate(tri, opts='p')
nnode = len(t['vertices'])
nelem = len(t['triangles'])
print("Génération d'un maillage avec {} noeuds et {} éléments".format(nnode, nelem))

print("~> Ecriture du maillage")

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

    # Ecriture du fichier XML
    with open(args.outfile_mesh, 'w') as fileout:
        fileout.write(template_render)

t2 = time.clock()
print("=> le temps d'execution est de : {}s".format(t2-t1))
