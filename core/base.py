"""
Classes :
* Coord
* ProfilTravers
* SuiteProfilsTravers
* LignesContrainte
* Lit
* MeshConstructor

Remarque : les abscisses recalculées par l'outil sont identiques à celle données par le module shapely (pour les LineString).
"""
from copy import deepcopy
from jinja2 import Environment, FileSystemLoader
from math import ceil
import numpy as np
from numpy.lib.recfunctions import append_fields
import os.path
import pandas as pd
from pyteltools.geom import BlueKenue as bk, Shapefile as shp
from pyteltools.geom import geometry
import shapefile
from shapely.geometry import LineString, Point
import sys
import time
import triangle

from .utils import get_intersections, float_vars, strictly_increasing


DIGITS = 4  # for csv and xml exports


class Coord:
    """
    Coord: ensemble de points 2D/3D consécutifs (ex: polylignes)

    ### Attributs
    - array

    ### Méthodes
    - compute_Xt
    - compute_xt
    - compute_xp
    - move_between_targets
    """
    def __init__(self, array, vars2add, remove_duplicates=False):
        """
        array: X, Y (Z is optional)
        Add Xt and xt columns if not already present
        """
        self.array = array
        vars = list(self.array.dtype.fields)
        self.z_labels = []
        for var in vars:
            if var not in ('X', 'Y', 'Xt', 'xt'):
                self.z_labels.append(var)

        if 'X' not in vars and 'Y' not in vars:
            sys.exit("Columns X and Y are compulsary. Following columns were found: {}".format(vars))

        if 'Xt' in vars2add:
            self.compute_Xt()
        if 'xt' in vars2add:
            self.compute_xt()

        if remove_duplicates:
            if not strictly_increasing(self.array['Xt']):
                print("ATTENTION : Des points doublons sont éliminés")
                # Suppression des doublons (points superposés dans la polyligne)
                points_a_conserver = np.ediff1d(self.array['Xt'], to_begin=1.) != 0.
                self.array = self.array[points_a_conserver]

    def compute_Xt(self):
        """
        Calcul de l'abscisse curviligne (distance 2D cumulée)
        """
        Xt = np.sqrt(np.power(np.ediff1d(self.array['X'], to_begin=0.), 2) +
                     np.power(np.ediff1d(self.array['Y'], to_begin=0.), 2))
        Xt = Xt.cumsum()
        self.array = append_fields(self.array, 'Xt', Xt, usemask=False)

    def compute_xt(self):
        """
        Calcul de l'asbcisse curviligne adimensionnée (de 0 à 1)
        /!\ La colonne 'Xt' doit déjà exister
        """
        if len(self.array) > 1:
            xt = (self.array['Xt'] - self.array['Xt'][0])/(self.array['Xt'][-1] - self.array['Xt'][0])
        else:
            xt = np.empty(len(self.array))
            xt.fill(-999.)
        self.array = append_fields(self.array, 'xt', xt, usemask=False)

    def move_between_targets(self, p1, p2):
        """
        @brief: Déplace la ligne pour la faire correspondre aux limites voulues
            (par translation pondérée linéairement de p1 à p2)
        @param p1 <Point>: point d'origine de la ligne
        @param p2 <Point>: point final de la ligne
        /!\ L'objet est modifié directement et la méthode ne retourne rien
        """
        array = deepcopy(self.array)
        xt = array['xt']
        self.array['X'] = array['X'] + (1-xt)*(p1.x-array['X'][0]) + xt*(p2.x-array['X'][-1])
        self.array['Y'] = array['Y'] + (1-xt)*(p1.y-array['Y'][0]) + xt*(p2.y-array['Y'][-1])

    def compute_xp(self):
        """Calcul de la distance projetée adimensionnée sur droite début->fin"""
        trace = LineString([self.array[['X', 'Y']][0], self.array[['X', 'Y']][-1]])
        Xp = np.empty(len(self.array), dtype=np.float)

        for i, row in enumerate(self.array):
            point = Point(list(row)[:2])
            Xp[i] = trace.project(point)

        if not strictly_increasing(Xp):
            sys.exit("L'abscisse projetée n'est pas strictement croissante le long du profil."
                     "Le profil n'est pas suffisamment droit...")

        xp = Xp/Xp[-1]
        self.array = append_fields(self.array, 'xp', xp, usemask=False)

    def convert_as_array(self):
        return np.column_stack((self.array['X'], self.array['Y'], self.array['Z']))

    def convert_as_linestring(self):
        return LineString(self.convert_as_array())

    def set_layers(self, z_values):
        self.z_labels = list(z_values.dtype.names)
        for z_label in self.z_labels:
            self.array = np.lib.recfunctions.append_fields(self.array, z_label, z_values[z_label], usemask=False)


class ProfilTravers:
    """
    ProfilTravers: représente un profil en travers

    ### Attributs
    - id <integer>: identifiant unique (numérotation automatique commençant à 0)
    - coord <Coord>: coordonnées (avec tous les niveaux)
    - geom <shapely.geometry.LineString>: objet géometrique 2D
    - limites <dict>: dictionnaire du type: {id_ligne: (Xt_profil, Xt_ligne, intersection.z)}
    - dist_proj_axe (créée par une méthode de <SuiteProfilsTravers>)

    ### Méthodes
    - _add_limit
    - get_limit_by_id
    - get_Xt_lignes
    - find_and_add_limit
    - common_limits
    - extraire_lit
    - sort_limites
    - calcul_nb_pts_inter
    - change_coord
    """
    LABELS = ['X', 'Y', 'Z', 'Xt']

    def __init__(self, id, coord, z_values, label):
        """
        Créer un profil à partir d'un semis de points ordonnés (coordonnées X,Y,Z)
        Aucune limite de lit n'est déclarée lorsque l'objet est créé
        """
        self.id = id
        self.label = label
        self.coord = Coord(np.array(coord, dtype=float_vars(['X', 'Y'])), ['Xt'])
        self.coord.set_layers(z_values)

        self.nb_points = len(self.coord.array)
        self.geom = LineString(coord)  # FIXME: contient toujours les points doublons

        self.limites = pd.DataFrame(columns=['Xt_profil', 'Xt_ligne', 'X', 'Y', 'Z', 'id_pt'])
        self.limites['id_pt'] = self.limites['id_pt'].astype(int)  # FIXME change dtype
        self.dist_proj_axe = -1

    def __repr__(self):
        return "{} #{} ({} points)".format(self.label, self.id, self.nb_points)

    def _add_limit(self, id_ligne, Xt_profil, Xt_ligne, point):
        """Ajoute une nouvelle limite au profil"""
        z_values = {}
        for label in self.coord.z_labels:
            z_values[label] = np.interp(Xt_profil, self.coord.array['Xt'], self.coord.array[label])
        row = pd.DataFrame({'Xt_profil': Xt_profil, 'Xt_ligne': Xt_ligne,
                            'X': point.x, 'Y': point.y, **z_values}, index=[id_ligne])
        self.limites = self.limites.append(row)

    def get_limit_by_id(self, id_ligne):
        return self.limites.loc[id_ligne]

    def get_Xt_lignes(self, id1, id2):
        return tuple(self.limites.loc[[id1, id2], 'Xt_ligne'])

    def find_and_add_limit(self, ligne_contrainte, dist_max=None):
        """
        @param ligne_contrainte <shapely.geometry.LineString>: ligne de contrainte 2D
        @param dist_max <float>: distance de tolérance pour détecter des intersections
        """
        if self.geom.intersects(ligne_contrainte.geom):
            intersection = self.geom.intersection(ligne_contrainte.geom)

            if isinstance(intersection, Point):
                # Calcul des projections
                Xt_profil = self.geom.project(intersection)
                Xt_ligne = ligne_contrainte.geom.project(intersection)
                self._add_limit(ligne_contrainte.id, Xt_profil, Xt_ligne, intersection)
            else:
                sys.exit("L'intersection entre '{}' et '{}' ne correspond pas à un point unique (type = {})".format(
                    self, ligne_contrainte, type(intersection)))

        else:
            distance = self.geom.distance(ligne_contrainte.geom)
            if dist_max is not None:
                if distance < dist_max:
                    # Test de proximité avec tous les points de la ligne de contrainte
                    for i, coord in enumerate(ligne_contrainte.coord):
                        point = Point(coord)
                        dist = self.geom.distance(point)
                        if dist < dist_max:
                            # On a trouvé le bon point et on ajoute la limite
                            Xt_ligne = ligne_contrainte.geom.project(point)
                            Xt_profil = self.geom.project(point)
                            intersection = self.geom.interpolate(Xt_profil)
                            self._add_limit(ligne_contrainte.id, Xt_profil, Xt_ligne, intersection)
                            print("Ajout de la limite avec la ligne {} après {} itérations (distance = {})".format(
                                ligne_contrainte.id, i, dist))
                            break

    def common_limits(self, other):
        """
        @brief: Liste les limites communes entre les deux profils (self et other)
        @param other <ProfilTravers>: profil en travers amont/aval
        """
        return self.limites.index[np.in1d(self.limites.index, other.limites.index, assume_unique=True)]

    def extraire_lit(self, id_lit_1, id_lit_2):
        """
        @brief: Extraire les coordonnées d'une partie du profil
            Les points du profil compris entre les deux lits et avec éventuellement les points de bord interpolés
            /!\ id_lit_1 et id_lit_2 doivent être dans l'ordre Xt croissant (plante sinon)
        @return <Lit>: avec les colonnes avec les colonnes ('X', 'Y', 'Z', 'Xt', 'xt')
        """
        limit1 = self.get_limit_by_id(id_lit_1)
        limit2 = self.get_limit_by_id(id_lit_2)

        Xt1 = limit1['Xt_profil']
        Xt2 = limit2['Xt_profil']

        # Vérifie que les Xt sont croissants du id_lit_1 au id_lit_2
        if Xt1 > Xt2:
            sys.exit("L'ordre des lits ({}, {}) n'est pas par Xt croissant pour le {}".format(id_lit_1, id_lit_2, self))

        Xt_profil = self.coord.array['Xt']
        sub_coord = self.coord.array[np.logical_and(Xt_profil >= Xt1, Xt_profil <= Xt2)]

        # Ajoute le premier point si nécessaire
        if Xt1 not in Xt_profil:
            row = np.array([tuple(limit1[var] if var != 'Xt' else Xt1 for var in sub_coord.dtype.names)],
                           dtype=sub_coord.dtype)
            sub_coord = np.insert(sub_coord, 0, row)

        # Ajoute le dernier point si nécessaire
        if Xt2 not in Xt_profil:
            row = np.array([tuple(limit2[var] if var != 'Xt' else Xt2 for var in sub_coord.dtype.names)],
                           dtype=sub_coord.dtype)
            sub_coord = np.append(sub_coord, row)

        # Vérification de l'ordre des points
        if not strictly_increasing(sub_coord['Xt']):
            print("/!\ Les Xt ne sont pas strictement croissants")  #FIXME: It should not append
            print(sub_coord['Xt'])
            print("Veuillez vérifier ci-dessus, avec les limites suivantes :")
            print(limit1)
            print(limit2)
            points_a_conserver = np.ediff1d(sub_coord['Xt'], to_begin=1.) != 0.
            sub_coord = sub_coord[points_a_conserver]

        return Lit(sub_coord, ['xt'])

    def sort_limites(self):
        """
        @brief: Trie les limites par Xt croissants
        Étape nécessaire pour chercher les limites communes entre deux profils plus facilement
        """
        self.limites = self.limites.sort_values(by='Xt_profil')

    def calcul_nb_pts_inter(self, other, pas_long):
        """
        @brief: Calcule le nombre de profils intermédiaires nécessaires
        /!\ Cette valeur englobe les points de bord (début et fin)
        @param other <ProfilTravers>: profil en travers
        @param pas_long <float>: pas longitudinal (en m)
        """
        dist_min_profil = self.geom.distance(other.geom)
        return ceil(dist_min_profil/pas_long) - 1

    def change_coord(self, array):
        self.coord.array = array
        self.coord.compute_Xt()
        self.nb_points = len(self.coord.array)
        self.geom = self.coord.convert_as_linestring()


def get_field_index(filename, field_id):
    if field_id is not None:
        names, _ = shp.get_attribute_names(filename)
        try:
            return names.index(field_id)
        except ValueError:
            sys.exit("Le champ `%s` n'existe pas" % field_id)


class SuiteProfilsTravers:
    """
    SuiteProfilsTravers: ensemble de profils en travers

    ### Attributs
    - suite <[ProfilTravers]>

    ### Méthodes
    - find_and_add_limits
    - calculer_dist_proj_axe
    """
    def __init__(self, filename, label, field_id=None):
        self.suite = []
        if filename.endswith('.i3s'):
            with bk.Read(filename) as in_i3s:
                in_i3s.read_header()
                for i, line in enumerate(in_i3s.get_open_polylines()):
                    line_id = i if field_id is None else line.attributes()[0]  # Use `Value` if field is not None
                    z_array = np.array([(coord[2], ) for coord in line.polyline().coords], dtype=float_vars('Z'))
                    if line.to_2d:
                        line = line.to_2d()
                    self.suite.append(ProfilTravers(line_id, list(line.polyline().coords), z_array, label))
        elif filename.endswith('.shp'):
            shp_type = shp.get_shape_type(filename)
            if shp_type in (shapefile.POLYLINEZ, shapefile.POLYLINEM):
                field_id_index = get_field_index(filename, field_id)
                for i, line in enumerate(shp.get_open_polylines(filename)):
                    line_id = i if field_id is None else line.attributes()[field_id_index]
                    z_array = np.array([(coord[2], ) for coord in line.polyline().coords], dtype=float_vars('Z'))
                    if line.to_2d:
                        line = line.to_2d()
                    self.suite.append(ProfilTravers(line_id, list(line.polyline().coords), z_array, label))
            elif shp_type == shapefile.POINTZ:
                field_id_index = get_field_index(filename, field_id)
                field_indexes, field_names = [], []
                for index, name in shp.get_numeric_attribute_names(filename):
                    if name.startswith('Z'):
                        field_indexes.append(index)
                        field_names.append(name)
                print('Variables : %s' % field_names)
                coords, z_layers = [], []
                last_point_id = None
                for i, (point, attributes) in enumerate(shp.get_points(filename, with_z=True)):
                    point_id = i if field_id is None else attributes[field_id_index]
                    if i > 0 and point_id != last_point_id:
                        z_array = np.array(z_layers, dtype=float_vars(['Z'] + field_names))
                        self.suite.append(ProfilTravers(last_point_id, coords, z_array, label))
                        coords, z_layers = [], []
                    coords.append(point[:2])
                    z_layers.append((point[2],) + tuple(attributes[index] for index in field_indexes))
                    last_point_id = point_id
                z_array = np.array(z_layers, dtype=float_vars(['Z'] + field_names))
                self.suite.append(ProfilTravers(last_point_id, coords, z_array, label))
            else:
                sys.exit("Le fichier %s n'est pas de type POINTZ ou POLYLINEZ[M]" % filename)

        else:
            raise NotImplementedError("Seuls les formats i3s et shp sont supportés pour les profils en travers")

    def __add__(self, other):
        newsuite = deepcopy(self)
        newsuite.suite = self.suite + other.suite
        return newsuite

    def __getitem__(self, index):
        if isinstance(index, slice):
            newsuite = deepcopy(self)
            newsuite.suite = self.suite[index]
            return newsuite
        else:
            return self.suite[index]

    def __len__(self):
        return len(self.suite)

    def __repr__(self):
        return [profil for profil in self.suite].__repr__()

    def find_and_add_limits(self, lignes_contraintes, dist_max):
        """
        @param lignes_contraintes <[LigneContrainte]>: lignes de contraintes
        @param dist_max <float>:
        """
        print("~> Recherche des limites de lits")
        for i, profil_travers in enumerate(self):
            for ligne_contrainte in lignes_contraintes:
                profil_travers.find_and_add_limit(ligne_contrainte, dist_max)
            profil_travers.sort_limites()
            limits = list(profil_travers.limites.index)  # only to print

            print("> {}".format(profil_travers))
            print("{} limites trouvées avec les lignes {}".format(len(limits), limits))

    def check_intersections(self):
        print("~> Vérifications non intersections des profils et des épis")
        intersections = get_intersections([profil.geom for profil in self])
        if intersections:
            print("Les intersections suivantes sont trouvées")
            for (i, j) in intersections:
                print("- entre '{}' et '{}'".format(self[i], self[j]))
            sys.exit("ERREUR: Des profils s'intersectent")

    def compute_dist_proj_axe(self, axe_geom):
        """
        @brief: Calculer la distance projetée sur l'axe
        @param axe_geom <shapely.geometry.LineString>: axe hydraulique
        /!\ Orientation de l'axe
        """
        print("~> Calcul des abscisses sur l'axe hydraulique (pour ordonner les profils/épis)")
        for profil in self:
            profil_geom = profil.geom
            if profil_geom.intersects(axe_geom):
                intersection = profil_geom.intersection(axe_geom)
                if isinstance(intersection, Point):
                    profil.dist_proj_axe = axe_geom.project(intersection)
                else:
                    sys.exit("L'intersection entre le '{}' et l'axe hydraulique n'est pas un point unique".format(profil))
            else:
                print(list(profil.geom.coords))
                sys.exit("ERREUR: Le '{}' n'intersection pas l'axe".format(profil))

    def sort_by_dist(self):
        self.suite = sorted(self.suite, key=lambda x: x.dist_proj_axe)

    def export_profil_shp(self, outfile_profils):
        """
        Write a shapefile with 3D Points
        @param outfile_profils <str>: output file name
        """
        w = shp.MyWriter(shapeType=shapefile.POINTZ)
        w.field('profil_id', 'C')
        w.field('Z', 'N', decimal=6)
        for profil in self.suite:
            for (x, y, z) in list(profil.geom.coords):
                w.point(x, y, z)
                w.record(str(profil.id), z)
        w.save(outfile_profils)


class LigneContrainte:
    """
    LigneContrainte: polyligne 2D ouverte permettant de distinguer les lits

    ### Attributs
    - id <integer>: identifiant unique (numérotation automatique commençant à 0)
    - coord <2D-array float>: coordonnées ['X', 'Y'] des points de la ligne
    - geom <shapely.geometry.LineString>: objet géométrique

    ### Méthodes
    - coord_sampling_along_line
    """
    def __init__(self, id, coord):
        """
        Création d'un objet LigneContrainte à partir des coordoonées de la ligne
        coord: liste de tuples (de taille 2)
          Ex: [(x1,y1), (x2,y2), ..., (xn,yn)]
        """
        self.id = id
        self.nb_points = len(coord)
        self.coord = np.array(coord)
        self.geom = LineString(coord)

    def __repr__(self):
        return "Ligne de contrainte n°{} ({} points)".format(self.id, self.nb_points)

    @staticmethod
    def get_lines_from_file(filename):
        """
        Extraction d'objects LineString à partir d'un fichier i2s
        Info: value is ignored
        """
        lines = []
        if filename is not None:
            if filename.endswith('.i2s'):
                with bk.Read(filename) as in_i2s:
                    in_i2s.read_header()
                    for i, line in enumerate(in_i2s.get_open_polylines()):
                        lines.append(LigneContrainte(i, list(line.polyline().coords)))
            elif filename.endswith('.shp'):
                if shp.get_shape_type(filename) not in (shapefile.POLYLINE, shapefile.POLYLINEZ, shapefile.POLYLINEM):
                    sys.exit("Le fichier %s n'est pas de type POLYLINE[ZM]" % filename)
                for i, line in enumerate(shp.get_open_polylines(filename)):
                    lines.append(LigneContrainte(i, list(line.polyline().coords)))
            else:
                raise NotImplementedError("Seuls les formats i2s et shp sont supportés pour les profils en travers")
        return lines

    @staticmethod
    def get_lines_from_profils(profils_en_travers):
        first_coords = []
        last_coords = []
        for profil in profils_en_travers:
            first_coords.append(profil.geom.coords[0][:2])
            last_coords.append(profil.geom.coords[-1][:2])
        return [LigneContrainte(0, first_coords), LigneContrainte(1, last_coords)]

    def coord_sampling_along_line(self, Xp1, Xp2, Xp_adm_int):
        """
        @brief: Extraire les coordonnées des points interpolés le long de la ligne entre les deux abscisses
        @param: Xp_adm_int <1D-array>: répartition des points entre Xp1 et Xp2 (valeurs de 0 à 1, avec un parcours de Xp1 à Xp2)
        @param: Xp1 <float>: abscisse de début
        @param: Xp2 <float>: abscisse de fin
        """
        # Construction de la liste des abscisses cibles (en mètres)
        Xp = (1 - Xp_adm_int)*Xp1 + Xp_adm_int*Xp2

        coord_int = []
        for Xp_cur in Xp:
            point = self.geom.interpolate(Xp_cur)
            coord_int.append(point.coords[0][:2])
        np_coord_int = np.array(coord_int, dtype=float_vars(['X', 'Y']))

        return np_coord_int


class Lit(Coord):
    """
    Lit(s): morceau de profil borné par des limites
    ~> Hérite de la classe <Coord>

    ### Attributs
    - coord
    - nb_pts_trans

    ### Méthodes
    Toutes ces méthodes retournent un array avec les colonnes VARS2INT
    - interp_along_lit_auto
    - interp_along_lit
    - interp_inter_lineaire
    """
    def var2int(self):
        return ['X', 'Y'] + self.z_labels

    def interp_along_lit_auto(self, pas_trans):
        nb_pts_trans = ceil((self.array['Xt'][-1] - self.array['Xt'][0])/pas_trans) + 1
        Xt_adm_list = np.linspace(0., 1., num=nb_pts_trans)
        return self.interp_along_lit(Xt_adm_list)

    def interp_along_lit(self, Xt_adm_list):
        """
        @brief: Echantionnage du lit (interpolation sur vecteur adimensionné de Xt)
        @param Xt_adm_list <1D-array float>: distance adimensionnée (valeurs entre 0 et 1)
        """
        nb_pts_trans = len(Xt_adm_list)
        array_ech = np.empty(nb_pts_trans, dtype=float_vars(self.var2int()))
        for label in self.var2int():
            array_ech[label] = np.interp(Xt_adm_list, self.array['xt'], self.array[label])
        return array_ech

    def interp_inter_lineaire(self, other, coeff, nb_pts_trans):
        """
        @brief: Interpolation linéaire entre deux lits
        @param other <ProfilTravers>: profil en travers amont ou aval
        @param coeff <float>: pondération entre self et other (0=self, 1=other)
        """
        Xt_adm_list = np.linspace(0., 1., num=nb_pts_trans)

        array_1 = self.interp_along_lit(Xt_adm_list)
        array_2 = other.interp_along_lit(Xt_adm_list)
        array = np.empty(nb_pts_trans, dtype=float_vars(self.var2int()))

        for var in self.var2int():
            array[var] = (1-coeff)*array_1[var] + coeff*array_2[var]
        return array


class MeshConstructor:
    """
    MeshConstructor: données pour construire un maillage

    ### Attributs
    - points
    - i_pt <int> (curseur pour répérer l'avancement)
    - segments

    ### Méthodes
    - add_points
    - add_segments
    - add_segments_from_node_list
    - export_as_dict

    /!\ Les éléments sont numérotés à partir de 0
      (dans le tableau segments)
    """
    POINTS_DTYPE = float_vars(['X', 'Y', 'profil']) + [(var, np.int) for var in ('lit', )]

    def __init__(self, profils_travers, pas_trans):
        self.profils_travers = profils_travers
        self.pas_trans = pas_trans

        self.z_labels = profils_travers[0].coord.z_labels
        self.points = np.empty(0, dtype=MeshConstructor.POINTS_DTYPE + float_vars(self.z_labels))
        self.i_pt = -1
        self.segments = np.empty([0, 2], dtype=np.int)
        self.triangle = None  # set by `build_mesh`

    def add_points(self, coord, profil, lit):
        """!
        @brief: Ajouter des sommets/noeuds
        @param coord <2D-array float>: tableau des coordonnées avec les colonnes ['X', 'Y', ...] (voir POINTS_DTYPE)
        """
        new_coord = np.empty(len(coord), dtype=self.points.dtype)
        for var in ['X', 'Y'] + self.z_labels:  # copy existing columns
            new_coord[var] = coord[var]
        #FIXME: avoid copying in using np.lib.recfunctions.append_fields?
        new_coord['profil'] = profil
        new_coord['lit'] = lit
        self.i_pt += len(new_coord)
        self.points = np.hstack((self.points, new_coord))

    def add_segments(self, seg):
        """!
        @brief: Ajouter des segments à contraintre
        @param seg <2D-array int>: série de couples de sommets
        """
        self.segments = np.hstack((self.segments, seg))

    def add_segments_from_node_list(self, node_list):
        """
        @brief: Ajoute les sommets à partir d'une liste de noeuds successifs
        @param node_list <1D-array int>: série de noeuds
        """
        new_segments = np.column_stack((node_list[:-1], node_list[1:]))
        self.segments = np.vstack((self.segments, new_segments))

    def export_as_dict(self):
        """
        @brief: Exporter les données pour triangle.triangulate
        """
        return {'vertices': np.array(np.column_stack((self.points['X'], self.points['Y']))),
                'segments': self.segments}

    def build_initial_profiles(self):
        print("~> Interpolation sur les profils existants en prenant en compte le passage des lignes de contraintes")
        for i in range(len(self.profils_travers)):
            cur_profil = self.profils_travers[i]
            print(cur_profil)

            # Recherche des limites communes "amont"
            if i == 0:
                common_limites_id_1 = cur_profil.limites.index
            else:
                common_limites_id_1 = cur_profil.common_limits(self.profils_travers[i - 1])

            # Recherche des limites communes "aval"
            if i == len(self.profils_travers) - 1:
                common_limites_id_2 = cur_profil.limites.index
            else:
                common_limites_id_2 = cur_profil.common_limits(self.profils_travers[i + 1])

            # Union (non ordonnée) des limites amont/aval
            limites_id = np.union1d(common_limites_id_1, common_limites_id_2)
            # Ré-ordonne les limites
            limites_id = cur_profil.limites.index[np.in1d(cur_profil.limites.index, limites_id, assume_unique=True)]

            first_lit = True
            for j, (id1, id2) in enumerate(zip(limites_id, limites_id[1:])):
                lit = cur_profil.extraire_lit(id1, id2)
                coord_int = lit.interp_along_lit_auto(self.pas_trans)

                if first_lit:
                    cur_profil.limites.loc[id1, 'id_pt'] = self.i_pt + 1
                else:
                    coord_int = coord_int[1:]

                self.add_points(coord_int, profil=cur_profil.dist_proj_axe, lit=j)

                cur_profil.limites.loc[id2, 'id_pt'] = self.i_pt

                # Ajoute les nouveaux segments
                new_i_pt = np.arange(cur_profil.limites['id_pt'].loc[id1],
                                     cur_profil.limites['id_pt'].loc[id2] + 1)
                self.add_segments_from_node_list(new_i_pt)

                if first_lit:
                    first_lit = False

    def build_interp(self, lignes_contraintes, pas_long, constant_ech_long):
        """
        @param profils_travers <SuiteProfilsTravers>:
        """
        self.build_initial_profiles()

        ### BOUCLE SUR L'ESPACE INTER-PROFIL
        print("~> Construction du maillage par zone interprofils puis par lit")
        for i, (prev_profil, next_profil) in enumerate(zip(self.profils_travers, self.profils_travers[1:])):
            print("> Zone n°{} : entre {} et {}".format(i, prev_profil, next_profil))

            if constant_ech_long:
                nb_pts_inter = prev_profil.calcul_nb_pts_inter(next_profil, pas_long)
                Xp_adm_list = np.linspace(0.0, 1.0, num=nb_pts_inter + 2)[1:-1]

            # Recherche des limites communes entre les deux profils
            common_limites_id = prev_profil.common_limits(next_profil)
            print("Limites de lits communes : {}".format(list(common_limites_id)))

            if len(common_limites_id) < 2:
                sys.exit("ERREUR: aucune interpolation pour l'intervalle {} ({} limites communes)".format(i, len(common_limites_id)))

            else:
                first_lit = True
                ### BOUCLE SUR LES LITS (= MORCEAU(X) DE PROFIL)
                for j, (id1, id2) in enumerate(zip(common_limites_id, common_limites_id[1:])):
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

                    if not constant_ech_long:
                        nb_pts_inter = ceil(min(dXp_L1, dXp_L2)/pas_long) - 1
                        Xp_adm_list = np.linspace(0.0, 1.0, num=nb_pts_inter + 2)[1:-1]

                    L1_coord_int = lignes_contraintes[id1].coord_sampling_along_line(Xp_profil1_L1, Xp_profil2_L1, Xp_adm_list)
                    L2_coord_int = lignes_contraintes[id2].coord_sampling_along_line(Xp_profil1_L2, Xp_profil2_L2, Xp_adm_list)

                    ### BOUCLE SUR LES LIGNES
                    for k in range(nb_pts_inter):
                        Xp = Xp_adm_list[k]
                        dist_proj_axe = prev_profil.dist_proj_axe*Xp + next_profil.dist_proj_axe*(1 - Xp)
                        P1 = Point(tuple(L1_coord_int[k]))
                        P2 = Point(tuple(L2_coord_int[k]))

                        lit_int = Lit(lit_1.interp_inter_lineaire(lit_2, Xp,
                                                                  ceil(P1.distance(P2)/self.pas_trans)+1), ['Xt', 'xt'])
                        lit_int.move_between_targets(P1, P2)
                        coord_int = lit_int.array[['X', 'Y'] + lit_int.z_labels]
                        pt_list_L1.append(self.i_pt + 1)

                        if not first_lit:
                            # ignore le 1er point car la ligne de contrainte a déjà été traitée
                            coord_int = coord_int[1:]

                        self.add_points(coord_int, profil=dist_proj_axe, lit=j)

                        pt_list_L2.append(self.i_pt)

                    pt_list_L2 = np.array([prev_profil.limites['id_pt'].loc[id2]] + pt_list_L2 + [next_profil.limites['id_pt'].loc[id2]])
                    self.add_segments_from_node_list(pt_list_L2)

                    if first_lit:
                        pt_list_L1 = np.array([prev_profil.limites['id_pt'].loc[id1]] + pt_list_L1 + [next_profil.limites['id_pt'].loc[id1]])
                        self.add_segments_from_node_list(pt_list_L1)
                        first_lit = False

    def corr_bathy_on_epis(self, epis, dist_corr_epi):
        print("~> Correction de la bathymétrie autour des épis")
        if set(self.z_labels) != set(['Z']):
            sys.exit('Impossible de corriger les épis sur les couches sédimentaires')
        for epi in epis:
            epi_geom = epi.coord.convert_as_linestring()
            for i, coord in enumerate(self.points):
                pt_node = Point(tuple(coord[:3]))
                if epi_geom.distance(pt_node) < dist_corr_epi:
                    Xt_proj = epi_geom.project(pt_node)
                    pt_proj = epi_geom.interpolate(Xt_proj)
                    self.points['Z'][i] = pt_proj.z

    def build_mesh(self):
        print("~> Calcul du maillage")
        tri = self.export_as_dict()
        self.triangle = triangle.triangulate(tri, opts='p')
        try:
            nnode, nelem = len(self.triangle['vertices']), len(self.triangle['triangles'])
        except KeyError:
            sys.exit("ERREUR: La génération du maillage a échouée!")
        print("Génération d'un maillage avec {} noeuds et {} éléments".format(nnode, nelem))

    def export_points(self, path):
        if path.endswith('.xyz'):
            print("~> Exports en xyz des points")
            with open(path, 'wb') as fileout:
                np.savetxt(fileout, self.points[['X', 'Y', 'Z']], fmt='%.4f')
        elif path.endswith('.shp'):
            pass
        else:
            raise NotImplementedError

    def export_profiles(self, path):
        """
        /!\ Pas cohérent si constant_ech_long est différent de True
        """
        # Export as points
        if path.endswith('_pts.shp'):
            w = shapefile.Writer(shapefile.POINTZ)
            # w.field('profil', 'C', '32')
            w.field('PK', 'N', decimal=6)
            # w.field('dist', 'N', decimal=6)
            for name in self.z_labels:
                w.field(name, 'N', decimal=6)
            for dist in np.unique(self.points['profil']):
                pos = self.points['profil'] == dist
                vars = self.points[pos].dtype.names
                for row in self.points[pos]:
                    values = {key: value for key, value in zip(vars, row)}
                    w.point(values['X'], values['Y'], values['Z'], shapeType=shapefile.POINTZ)
                    w.record(dist, *[values[name] for name in self.z_labels])
            w.save(path)
            return

        # Export as lines
        lines = []
        for dist in np.unique(self.points['profil']):
            pos = self.points['profil'] == dist
            line = geometry.Polyline([(x, y, z) for x, y, z in self.points[pos][['X', 'Y', 'Z']]])
            line.add_attribute(dist)
            lines.append(line)

        if path.endswith('.i3s'):
            with bk.Write(path) as out_i3s:
                out_i3s.write_header()
                out_i3s.write_lines(lines, line.attributes())
        elif path.endswith('.shp'):
            shp.write_shp_lines(path, shapefile.POLYLINEZ, lines, 'Z')
        else:
            raise NotImplementedError

    def export_mesh(self, path):
        print("~> Écriture du maillage")

        nnode, nelem = len(self.triangle['vertices']), len(self.triangle['triangles'])
        if path.endswith(".t3s"):
            with open(path, 'w', newline='') as fileout:
                # Écriture en-tête
                date = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
                fileout.write("""#########################################################################
:FileType t3s  ASCII  EnSim 1.0
# Canadian Hydraulics Centre/National Research Council (c) 1998-2012
# DataType                 2D T3 Scalar Mesh
#
:Application              BlueKenue
:Version                  3.3.4
:WrittenBy                mailleurtatooine
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

            with open(path, mode='ab') as fileout:
                # Tableau des coordonnées (x, y, z)
                np.savetxt(fileout, np.column_stack((self.triangle['vertices'], self.points['Z'])), delimiter=' ',
                           fmt='%.{}f'.format(DIGITS))

                # Tableau des éléments (connectivité)
                np.savetxt(fileout, self.triangle['triangles'] + 1, delimiter=' ', fmt='%i')

        elif path.endswith(".xml"):
            env = Environment(
                loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'data')))
            template = env.get_template("LandXML_template.xml")
            template_render = template.render(
                nodes=np.round(np.column_stack((self.triangle['vertices'], self.points['Z'])), DIGITS),
                ikle=self.triangle['triangles'] + 1
            )

            # Écriture du fichier XML
            with open(path, 'w') as fileout:
                fileout.write(template_render)

        else:
            raise NotImplementedError
