"""
Classes :
* Coord
* ProfilTravers
* SuiteProfilsTravers
* LigneContrainte
* Lit
* MeshConstructor

Remarque : les abscisses recalculées par l'outil sont identiques à celle données par le module shapely (pour les LineString).
"""
from collections import OrderedDict
from copy import deepcopy
from jinja2 import Environment, FileSystemLoader
import math
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib.recfunctions import append_fields, rename_fields
import os.path
from pyteltools.geom import BlueKenue as bk, Shapefile as shp
from pyteltools.geom import geometry
from pyteltools.slf import Serafin
from pyteltools.slf.variable.variables_2d import basic_2D_vars_IDs
from scipy import interpolate
import shapefile
from shapely.geometry import LineString, MultiPoint, Point
import time
import triangle

from .spline_cubic_hermite import CubicHermiteSpline
from .utils import float_vars, get_intersections, logger, strictly_increasing, TatooineException



DIGITS = 4  # for csv and xml exports
LANG = 'fr'  # for variable names
COURLIS_FLOAT_FMT = '%.6f'


class Coord:
    """
    Coord: ensemble de points 2D/3D consécutifs (ex: polylignes)

    ### Attributs
    - array: structured array with coordinates ('X', 'Y'), variables (see z_labels) and eventually distance(s) ('Xt', 'xt')
    - z_labels: nom des variables d'intérêt
    - values: structured array with values (different variables are possible)

    ### Méthodes
    - compute_Xt
    - compute_xt
    - move_between_targets
    - compute_xp
    - convert_as_array
    - convert_as_linestring
    """
    XY = ['X', 'Y']

    def __init__(self, array, vars2add, remove_duplicates=False):
        """
        @param array: structured array. `X`, `Y` are compulsory but other variables (such as Z) and distances are optional
        @param var2add <[str]>: list including eventually `Xt` and/or `xt` top compute them if they are not already present
        @param remove_duplicates <bool>: remove consecutive duplicated points
        """
        self.array = array
        self.coord_labels = list(self.array.dtype.fields.keys())
        self.values = None

        if 'X' not in self.coord_labels and 'Y' not in self.coord_labels:
            raise TatooineException("Columns X and Y are compulsory.")

        if 'Xt' in vars2add:
            self.compute_Xt()
        if 'xt' in vars2add:
            self.compute_xt()

        if remove_duplicates:
            if not strictly_increasing(self.array['Xt']):
                logger.warn("Des points doublons sont éliminés")
                # Suppression des doublons (points superposés dans la polyligne)
                points_a_conserver = np.ediff1d(self.array['Xt'], to_begin=1.) != 0.
                self.array = self.array[points_a_conserver]

    def nb_var(self):
        if self.values is None:
            return 0
        else:
            return len(self.values.dtype.names)

    def compute_Xt(self):
        """
        Calcul de l'abscisse curviligne `Xt` (distance 2D cumulée)
        """
        Xt = np.sqrt(np.power(np.ediff1d(self.array['X'], to_begin=0.), 2) +
                     np.power(np.ediff1d(self.array['Y'], to_begin=0.), 2))
        Xt = Xt.cumsum()

        # Update or append `Xt` array column
        if 'Xt' in self.array.dtype.fields:
            self.array['Xt'] = Xt
        else:
            self.array = append_fields(self.array, 'Xt', Xt, usemask=False)

    def compute_xt(self):
        """
        Calcul de l'abscisse curviligne adimensionnée `xt` (de 0 à 1)
        /!\ La colonne `Xt` doit déjà exister
        """
        if len(self.array) > 1:
            xt = (self.array['Xt'] - self.array['Xt'][0])/(self.array['Xt'][-1] - self.array['Xt'][0])
        else:
            xt = np.empty(len(self.array))
            xt.fill(-999.)

        # Update or append `xt` array column
        if 'xt' in self.array.dtype.fields:
            self.array['xt'] = xt
        else:
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
            raise TatooineException("L'abscisse projetée n'est pas strictement croissante le long du profil. "
                                    "Le profil n'est pas suffisamment droit...")

        xp = Xp/Xp[-1]
        self.array = append_fields(self.array, 'xp', xp, usemask=False)

    def convert_as_array(self):
        """Returns a float 2D-array with X, Y and Z"""
        return np.column_stack((self.array['X'], self.array['Y'], self.values['Z']))

    def convert_as_linestring(self):
        """Returns a LineString object"""
        return LineString(self.convert_as_array())


class ProfilTravers:
    """
    ProfilTravers: représente un profil en travers avec éventuellement plusieurs variables en chacun de ses points

    ### Attributs
    - id <integer|str>: identifiant unique
    - label <str>: type de profil (`Profils en travers` ou `Épi`)
    - coord <Coord>: coordonnées (X et Y) avec toutes les variables
    - nb_points <int>: number of points
    - geom <shapely.geometry.LineString>: objet géometrique 2D
    - limites <OrderedDict>: dictionnaire ordonné avec comme clé `id_ligne` et comme valeur
          un dict avec Xt_profil Xt_ligne, X, Y, ...
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
    def __init__(self, id, coord, label='Profil en travers'):
        """
        Créer un profil à partir d'un semis de points ordonnés (coordonnées X, Y)
        Aucune limite de lit n'est déclarée lorsque l'objet est créé

        @param id <integer|str>: identifiant unique
        @param label <str>: type de profil (`Profil en travers` ou `Épi`)
        @param coords <[tuple]>: sequence of X and Y coordinates
        """
        self.id = id
        self.label = label
        self.coord = Coord(np.array(coord, dtype=float_vars(['X', 'Y'])), ['Xt', 'xt'])
        self.nb_points = len(self.coord.array)
        self.geom = LineString(coord)  # FIXME: contient toujours les points doublons

        self.limites = OrderedDict()
        self.dist_proj_axe = -1

    def __repr__(self):
        return "{} #{} ({} points)".format(self.label, self.id, self.nb_points)

    def _add_limit(self, id_ligne, Xt_profil, Xt_ligne, point):
        """Ajoute une nouvelle limite au profil"""
        z_values = {}
        for label in self.coord.coord_labels:
            z_values[label] = np.interp(Xt_profil, self.coord.array['Xt'], self.coord.array[label])
        self.limites[id_ligne] = {'Xt_profil': Xt_profil, 'Xt_ligne': Xt_ligne,
                                  'X': point.x, 'Y': point.y, **z_values}

    def get_limit_by_id(self, id_ligne):
        return self.limites[id_ligne]

    def get_limit_by_idx(self, idx):
        return self.limites[list(self.limites.keys())[idx]]

    def get_Xt_lignes(self, id1, id2):
        return self.get_limit_by_id(id1)['Xt_ligne'], self.get_limit_by_id(id2)['Xt_ligne']

    def find_and_add_limit(self, ligne_contrainte, dist_max=None):
        """
        @param ligne_contrainte <shapely.geometry.LineString>: ligne de contrainte 2D
        @param dist_max <float>: distance de tolérance pour détecter des intersections
        """
        if self.geom.intersects(ligne_contrainte.geom):
            intersection = self.geom.intersection(ligne_contrainte.geom)

            if isinstance(intersection, MultiPoint):
                logger.warn("L'intersection entre '{}' et '{}' contient plusieurs points, "
                            "seul le premier est gardé.".format(self, ligne_contrainte))
                intersection = intersection[0]
            if isinstance(intersection, Point):
                # Calcul des projections
                Xt_profil = self.geom.project(intersection)
                Xt_ligne = ligne_contrainte.geom.project(intersection)
                self._add_limit(ligne_contrainte.id, Xt_profil, Xt_ligne, intersection)

            else:
                raise TatooineException("L'intersection entre '{}' et '{}' ne correspond pas rien!".format(
                    self, ligne_contrainte))
        else:
            if dist_max is not None:
                distance = self.geom.distance(ligne_contrainte.geom)
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
                            logger.debug("Ajout de la limite avec la ligne {} après {} itérations (distance = {})"
                                         .format(ligne_contrainte.id, i, dist))
                            break

    def common_limits(self, limite_ids):
        """
        @brief: Liste les limites communes entre les deux profils (self et other)
        @param limites
        """
        out_limites = []
        for limite_id in self.limites.keys():
            if limite_id in limite_ids:
                out_limites.append(limite_id)
        return out_limites

    def extraire_lit(self, id_lit_1, id_lit_2):
        """
        @brief: Extraire les coordonnées d'une partie du profil
            Les points du profil compris entre les deux lits et avec éventuellement les points de bord interpolés
            /!\ id_lit_1 et id_lit_2 doivent être dans l'ordre Xt croissant (plante sinon)
        @return <Lit>: avec les colonnes avec les colonnes ('X', 'Y', 'Xt', 'xt')
        """
        limit1 = self.get_limit_by_id(id_lit_1)
        limit2 = self.get_limit_by_id(id_lit_2)

        Xt1 = limit1['Xt_profil']
        Xt2 = limit2['Xt_profil']

        # Vérifie que les Xt sont croissants du id_lit_1 au id_lit_2
        if Xt1 > Xt2:
            raise TatooineException("L'ordre des lits ({}, {}) n'est pas par Xt croissant pour le {}".format(
                id_lit_1, id_lit_2, self))

        Xt_profil = self.coord.array['Xt']
        sub_coord = self.coord.array[np.logical_and(Xt_profil >= Xt1, Xt_profil <= Xt2)]

        # Ajoute le premier point si nécessaire
        if Xt1 not in Xt_profil:
            row = np.array([tuple(limit1[var] if var not in ('Xt', 'xt') else Xt1 for var in sub_coord.dtype.names)],
                           dtype=sub_coord.dtype)
            row['xt'] = 0.0
            sub_coord = np.insert(sub_coord, 0, row)

        # Ajoute le dernier point si nécessaire
        if Xt2 not in Xt_profil:
            row = np.array([tuple(limit2[var] if var not in ('Xt', 'xt') else Xt2 for var in sub_coord.dtype.names)],
                           dtype=sub_coord.dtype)
            row['xt'] = 1.0
            sub_coord = np.append(sub_coord, row)

        # Vérification de l'ordre des points
        if not strictly_increasing(sub_coord['Xt']):
            logger.debug("/!\ Les Xt ne sont pas strictement croissants")  # FIXME: It should not append
            logger.debug(sub_coord['Xt'])
            logger.debug("Veuillez vérifier ci-dessus, avec les limites suivantes :")
            logger.debug(limit1)
            logger.debug(limit2)
            points_to_keep = np.ediff1d(sub_coord['Xt'], to_begin=1.) != 0.
            sub_coord = sub_coord[points_to_keep]

        return Lit(sub_coord, ['xt'])

    def sort_limites(self):
        """
        @brief: Trie les limites par Xt croissants
        Étape nécessaire pour chercher les limites communes entre deux profils plus facilement
        """
        self.limites = OrderedDict(sorted(self.limites.items(), key=lambda x: x[1]['Xt_profil']))

    def calcul_nb_pts_inter(self, other, pas_long):
        """
        @brief: Calcule le nombre de profils intermédiaires nécessaires
        /!\ Cette valeur englobe les points de bord (début et fin)
        @param other <ProfilTravers>: profil en travers
        @param pas_long <float>: pas longitudinal (en m)
        """
        dist_min_profil = self.geom.distance(other.geom)
        return math.ceil(dist_min_profil/pas_long) - 1

    def project_straight_line(self):
        """Planar projection on a straight line joining first and ending point of the cross-section"""
        if self.limites:
            raise TatooineException("Limits have to be set after calling to project_straight_line!")

        # Build straight line
        first_point = (self.coord.array['X'][0], self.coord.array['Y'][0])
        last_point = (self.coord.array['X'][-1], self.coord.array['Y'][-1])
        line = LineString((first_point, last_point))

        # Update X, Y, Xt columns in self.coord.array
        for row in self.coord.array:
            point = Point((row['X'], row['Y']))
            dist = line.project(point)
            point_project = line.interpolate(dist)
            row['X'] = point_project.x
            row['Y'] = point_project.y
        self.coord.compute_Xt()

        # Update geom
        self.geom = self.coord.convert_as_linestring()

    def get_segments(self):
        x = self.coord.array['Xt']
        z = self.coord.values['Z']
        return [(xx, zz) for xx, zz in zip(x, z)]

    def get_angles(self):
        """Angles are within range [0, 360]"""
        angles = []
        for i, ((x1, z1), (x2, z2), (x3, z3)) in enumerate(zip(self.get_segments(), self.get_segments()[1:],
                                                               self.get_segments()[2:])):
            angle = math.degrees(math.atan2(z3 - z2, x3 - x2) - math.atan2(z1 - z2, x1 - x2))
            if angle < 0:
                angle += 360
            angles.append(angle)
        return angles

    def export_plot_travers(self, fig_path, overwrite=False):
        x = self.coord.array['Xt']
        z = self.coord.values['Z']

        fig, ax1 = plt.subplots(figsize=(16, 9))

        ax1.set_xlabel('Xt (m)')

        color = 'tab:red'
        ax1.set_ylabel('Z', color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.plot(x, z, marker='o', color=color, label='Z')

        color = 'tab:green'
        for limit_name, limit in self.limites.items():
            ax1.axvline(x=limit['Xt_profil'], linestyle='-', color=color)

        color = 'tab:blue'
        ax2 = ax1.twinx()
        ax2.set_ylabel('Angles (°)', color=color)
        ax2.plot(x[1:-1], np.array(self.get_angles()) - 180.0, color=color, label='Angles')
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.set_ylim(-180, 180)
        plt.yscale('symlog')

        if overwrite and os.path.exists(fig_path):
            os.remove(fig_path)
        plt.legend(loc='upper center')
        plt.savefig(fig_path, dpi=400)


def get_field_index(filename, field_id):
    if field_id is not None:
        names, _ = shp.get_attribute_names(filename)
        try:
            return names.index(field_id)
        except ValueError:
            raise TatooineException("Le champ `%s` n'existe pas" % field_id)


class SuiteProfilsTravers:
    """
    SuiteProfilsTravers: ensemble de profils en travers

    ### Attributs
    - suite <[ProfilTravers]>

    ### Méthodes
    - find_and_add_limits
    - calculer_dist_proj_axe
    """

    def __init__(self):
        self.suite = []

    def add_profile(self, profile):
        self.suite.append(profile)

    @staticmethod
    def from_file(filename, label, field_id=None, project_straight_line=False):
        profils_travers = SuiteProfilsTravers()
        if filename.endswith('.i3s'):
            with bk.Read(filename) as in_i3s:
                in_i3s.read_header()
                for i, line in enumerate(in_i3s.get_open_polylines()):
                    line_id = i if field_id is None else line.attributes()[0]  # Use `Value` if field is not None
                    z_array = np.array([(coord[2],) for coord in line.polyline().coords], dtype=float_vars('Z'))
                    line = line.to_2d()
                    profile = ProfilTravers(line_id, list(line.polyline().coords), label)
                    profile.coord.values = z_array
                    profils_travers.add_profile(profile)
        elif filename.endswith('.shp'):
            shp_type = shp.get_shape_type(filename)
            if shp_type in (shapefile.POLYLINEZ, shapefile.POLYLINEM):
                field_id_index = get_field_index(filename, field_id)
                for i, line in enumerate(shp.get_open_polylines(filename)):
                    line_id = i if field_id is None else line.attributes()[field_id_index]
                    z_array = np.array([(coord[2],) for coord in line.polyline().coords], dtype=float_vars(['Z']))
                    line = line.to_2d()
                    profile = ProfilTravers(line_id, list(line.polyline().coords), label)
                    profile.coord.values = z_array
                    profils_travers.add_profile(profile)
            elif shp_type == shapefile.POINTZ:
                field_id_index = get_field_index(filename, field_id)
                field_indexes, field_names = [], []
                for index, name in shp.get_numeric_attribute_names(filename):
                    if name.startswith('Z'):
                        field_indexes.append(index)
                        field_names.append(name)
                coords, z_layers = [], []
                last_point_id = None
                for i, (point, attributes) in enumerate(shp.get_points(filename, with_z=True)):
                    point_id = attributes[field_id_index]  # FIXME: should raise exception if field_id_index is None!
                    if i > 0 and point_id != last_point_id:
                        z_array = np.array(z_layers, dtype=float_vars(['Z'] + field_names))
                        profile = ProfilTravers(last_point_id, coords, label)
                        profile.coord.values = z_array
                        profils_travers.add_profile(profile)
                        coords, z_layers = [], []
                    coords.append(point[:2])
                    z_layers.append((point[2],) + tuple(attributes[index] for index in field_indexes))
                    last_point_id = point_id
                z_array = np.array(z_layers, dtype=float_vars(['Z'] + field_names))
                profile = ProfilTravers(last_point_id, coords, label)
                profile.coord.values = z_array
                profils_travers.add_profile(profile)
            else:
                raise TatooineException("Le fichier %s n'est pas de type POINTZ ou POLYLINEZ[M]" % filename)
        else:
            raise NotImplementedError("Seuls les formats i3s, shp sont supportés pour les profils en travers")

        if project_straight_line:
            for profile in profils_travers:
                profile.project_straight_line()
        return profils_travers

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
        @param dist_max <float>: distance de tolérance pour détecter des intersections
        """
        logger.info("~> Recherche des limites de lits")
        for i, profil_travers in enumerate(self):
            for ligne_contrainte in lignes_contraintes:
                profil_travers.find_and_add_limit(ligne_contrainte, dist_max)
            profil_travers.sort_limites()
            limits = profil_travers.limites.keys()  # only to print

            logger.debug("> {}".format(profil_travers))
            logger.debug("{} limites trouvées avec les lignes {}".format(len(limits), list(limits)))

    def check_intersections(self):
        """Display warning with intersections details"""
        logger.info("~> Vérifications non intersections des profils et des épis")
        intersections = get_intersections([profil.geom for profil in self])
        if intersections:
            logger.warn("Les intersections suivantes sont trouvées")
            for (i, j) in intersections:
                logger.warn("- entre '{}' et '{}'".format(self[i], self[j]))

    def compute_dist_proj_axe(self, axe_geom, dist_max):
        """
        @brief: Calculer la distance projetée sur l'axe
        @param axe_geom <shapely.geometry.LineString>: axe hydraulique (/!\ Orientation de l'axe)
        @param dist_max <float>: distance de tolérance pour détecter des intersections
        """
        logger.info("~> Calcul des abscisses sur l'axe hydraulique (pour ordonner les profils/épis)")
        to_keep_list = []
        for profil in self:
            profil_geom = profil.geom
            if profil_geom.intersects(axe_geom):
                intersection = profil_geom.intersection(axe_geom)
                if isinstance(intersection, Point):
                    profil.dist_proj_axe = axe_geom.project(intersection)
                else:
                    raise TatooineException("L'intersection entre le '{}' et l'axe hydraulique "
                                            "n'est pas un point unique".format(profil))
            else:
                if dist_max is not None:
                    for pos in (0, -1):
                        dist = profil_geom.distance(Point(axe_geom.coords[pos]))
                        if dist < dist_max:
                            profil.dist_proj_axe = 0.0 if pos == 0 else axe_geom.length
            if profil.dist_proj_axe == -1:
                logger.warn("{} n'intersecte pas l'axe (distance = {}m) et est ignoré".format(
                    profil, profil.geom.distance(axe_geom)))
                to_keep_list.append(False)
            else:
                to_keep_list.append(True)

        self.suite = [p for p, to_keep in zip(self.suite, to_keep_list) if to_keep]

    def sort_by_dist(self):
        self.suite = sorted(self.suite, key=lambda x: x.dist_proj_axe)

    def export_profil_shp(self, outfile_profils):
        """
        Write a shapefile with 3D LineString
        @param outfile_profils <str>: output file name
        """
        with shapefile.Writer(outfile_profils, shapeType=shapefile.POLYLINEZ) as w:
            w.field('profil_id', 'C')
            for profil in self.suite:
                array = profil.coord.array
                z_array = profil.coord.values['Z']
                coords = [(row['X'], row['Y'], z) for row, z in zip(array, z_array)]
                w.linez([coords])
                w.record(profil_id=str(profil.id))

    def add_constant_layer(self, name, thickness):
        for profil in self.suite:
            coord = profil.coord
            coord.add_single_layer(name, coord.array['Z'] - thickness)


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
    def __init__(self, id, coord, interp_coord='LINEAR'):
        """
        Création d'un objet LigneContrainte à partir des coordoonées de la ligne
        coord: liste de tuples (de taille 2)
          Ex: [(x1, y1), (x2, y2), ..., (xn, yn)]
        """
        self.id = id
        self.nb_points = len(coord)
        self.coord = np.array(coord)
        Xt = np.sqrt(np.power(np.ediff1d(self.coord[:, 0], to_begin=0.), 2) +
                     np.power(np.ediff1d(self.coord[:, 1], to_begin=0.), 2))
        self.Xt = Xt.cumsum()
        self.geom = LineString(coord)

        if interp_coord == 'LINEAR':
            self.interp = self.build_interp_linear()
        else:
            if interp_coord == 'CARDINAL':
                tan_method = CubicHermiteSpline.CARDINAL
            elif interp_coord == 'FINITE_DIFF':
                tan_method = CubicHermiteSpline.FINITE_DIFF
            else:
                raise NotImplementedError
            self.interp = self.build_interp_chs(tan_method)

    def __repr__(self):
        return "Ligne de contrainte n°{} ({} points)".format(self.id, self.nb_points)

    @staticmethod
    def get_lines_from_file(filename, interp_coord='LINEAR'):
        """
        Extraction d'objects LineString à partir d'un fichier i2s
        Info: value is ignored
        TODO 2: Use a SuiteProfilTravers instance to get the correct profiles order?
        """
        lines = []
        if filename is not None:
            if filename.endswith('.i2s'):
                with bk.Read(filename) as in_i2s:
                    in_i2s.read_header()
                    for i, line in enumerate(in_i2s.get_open_polylines()):
                        lines.append(LigneContrainte(i, list(line.polyline().coords), interp_coord))
            elif filename.endswith('.shp'):
                if shp.get_shape_type(filename) not in (shapefile.POLYLINE, shapefile.POLYLINEZ, shapefile.POLYLINEM):
                    raise TatooineException("Le fichier %s n'est pas de type POLYLINE[ZM]" % filename)
                for i, line in enumerate(shp.get_open_polylines(filename)):
                    lines.append(LigneContrainte(i, list(line.polyline().coords), interp_coord))
            else:
                raise NotImplementedError("Seuls les formats i2s et shp sont supportés pour les lignes de contrainte")
        return lines

    @staticmethod
    def get_lines_and_set_limits_from_profils(profils_en_travers, interp_coord='LINEAR'):
        # Build lines from cross-section extremities
        first_coords = []
        last_coords = []
        for profile in profils_en_travers:
            first_coords.append(profile.geom.coords[0][:2])
            last_coords.append(profile.geom.coords[-1][:2])
        lines = [LigneContrainte(0, first_coords, interp_coord),
                 LigneContrainte(1, last_coords, interp_coord)]

        # Set limits
        for id_line, line in enumerate(lines):
            for profile, Xt_ligne in zip(profils_en_travers, line.Xt):
                Xt_profil = profile.coord.array['Xt'][0] if id_line == 0 else profile.coord.array['Xt'][-1]
                intersection = profile.geom.interpolate(Xt_profil)
                profile._add_limit(id_line, Xt_profil, Xt_ligne, intersection)

        return lines

    def build_interp_linear(self):
        def interp_xy_linear(Xt_new):
            coord_int = []
            for dist in Xt_new:
                point = self.geom.interpolate(dist)
                coord_int.append(point.coords[0][:2])
            return np.array(coord_int, dtype=float_vars(['X', 'Y']))
        return interp_xy_linear

    def build_interp_chs(self, tan_method):
        spline_x = CubicHermiteSpline()  # x = spline_x(Xt)
        spline_y = CubicHermiteSpline()  # y = spline_y(Xt)

        spline_x.Initialize(np.vstack((self.Xt, self.coord[:, 0])).T, tan_method=tan_method)
        spline_y.Initialize(np.vstack((self.Xt, self.coord[:, 1])).T, tan_method=tan_method)

        def interp_xy_chs(Xt_new):
            coord_int = []
            for dist in Xt_new:
                x = spline_x.evaluate(dist)
                y = spline_y.evaluate(dist)
                coord_int.append((x, y))
            np_coord_int = np.array(coord_int, dtype=float_vars(['X', 'Y']))
            return np_coord_int
        return interp_xy_chs

    def coord_sampling_along_line(self, Xp1, Xp2, Xp_adm_int):
        """
        @brief: Extraire les coordonnées des points interpolés le long de la ligne entre les deux abscisses
        @param: Xp_adm_int <1D-array>: répartition des points entre Xp1 et Xp2 (valeurs de 0 à 1, avec un parcours de Xp1 à Xp2)
        @param: Xp1 <float>: abscisse de début
        @param: Xp2 <float>: abscisse de fin
        """
        # Construction de la liste des abscisses cibles (en mètres)
        Xp = (1 - Xp_adm_int)*Xp1 + Xp_adm_int*Xp2
        return self.interp(Xp)


class Lit(Coord):
    """
    Lit(s): morceau de profil borné par des limites
    ~> Hérite de la classe <Coord>

    ### Attributs
    - coord

    ### Méthodes
    Toutes ces méthodes retournent un array avec les colonnes VARS2INT
    - interp_along_lit_auto
    - interp_along_lit
    - interp_coord_linear
    """
    def interp_coord_along_lit_auto(self, pas_trans, nb_pts_trans=None):
        if nb_pts_trans is None:
            nb_pts_trans = math.ceil((self.array['Xt'][-1] - self.array['Xt'][0])/pas_trans) + 1
        Xt_adm_list = np.linspace(0., 1., num=nb_pts_trans)
        return self.interp_coord_along_lit(Xt_adm_list)

    def interp_coord_along_lit(self, Xt_adm_list):
        """
        @brief: Echantionnage du lit (interpolation sur vecteur adimensionné de Xt)
        @param Xt_adm_list <1D-array float>: distance adimensionnée (valeurs entre 0 et 1)
        """
        nb_pts_trans = len(Xt_adm_list)
        array = np.empty(nb_pts_trans, dtype=float_vars(Coord.XY + ['Xt', 'xt']))
        for label in Coord.XY + ['Xt']:
            array[label] = np.interp(Xt_adm_list, self.array['xt'], self.array[label])
        array['xt'] = Xt_adm_list
        return array

    def interp_coord_linear(self, other, coeff, nb_pts_trans):
        """
        @brief: Interpolation linéaire entre deux lits
        @param other <Lit>: profil en travers amont ou aval
        @param coeff <float>: pondération entre self et other (0=self, 1=other)
        @param nb_pts_trans <int>: nombre de points transversalement
        """
        Xt_adm_list = np.linspace(0., 1., num=nb_pts_trans)

        array_1 = self.interp_coord_along_lit(Xt_adm_list)
        array_2 = other.interp_coord_along_lit(Xt_adm_list)

        array = np.empty(nb_pts_trans, dtype=float_vars(Coord.XY + ['xt', 'Xt_amont', 'Xt_aval']))
        for var in Coord.XY:
            array[var] = (1 - coeff)*array_1[var] + coeff*array_2[var]
        array['Xt_amont'] = array_1['Xt']
        array['Xt_aval'] = array_2['Xt']
        return array


class MeshConstructor:
    """
    MeshConstructor: construire un maillage et interpoler les valeurs

    ### Attributs
    - profils_travers <SuiteProfilsTravers>
    - pas_trans <float>: pas transversal (m)
    - nb_pts_trans <int>: nombre de noeuds transversalement
    #TODO
    - points: structured array with columns ['X', 'Y', 'Xt_amont', 'Xt_aval', 'Xl', 'xl', 'zone', 'lit']
    - i_pt <int>: curseur pour repérer l'avancement
    - segments <2D-array int>: list of nodes numbers (0-indexed) to define constrainted segments
    - triangle <dict>: dictionary with 2 keys `vertices` (2D nodes) and `triangles` (connectivity table)

    ### Méthodes
    - var_names
    - add_points
    - add_segments
    - add_segments_from_node_list
    - export_triangulation_dict
    - export_floworiented_triangulation_dict
    - build_initial_profiles
    - build_interp
    - corr_bathy_on_epis
    - build_mesh
    - summary
    - export_points
    - export_segments
    - export_profiles
    - export_mesh
    - interp_values_from_profiles
    - get_merge_triangulation
    """
    POINTS_DTYPE = float_vars(['X', 'Y', 'xt', 'Xt_amont', 'Xt_aval', 'Xl', 'xl']) + \
                              [(var, np.int) for var in ('zone', 'lit')]

    def __init__(self, profils_travers, pas_trans=None, nb_pts_trans=None, interp_values='LINEAR'):
        self.profils_travers = profils_travers
        self.pas_trans = pas_trans
        self.nb_pts_trans = nb_pts_trans
        self.interp_values = interp_values

        self.points = np.empty(0, dtype=MeshConstructor.POINTS_DTYPE)
        self.nodes_values = np.empty([0, 0, 0], dtype=np.float)
        self.i_pt = int(-1)
        self.segments = np.empty([0, 2], dtype=np.int)
        self.triangle = None  # set by `build_mesh`

    def var_names(self):
        return list(self.profils_travers[0].coord.values.dtype.names)

    def add_points(self, coord, xl, index_zone, index_lit):
        """!
        @brief: Ajouter des sommets/noeuds
        @param coord <2D-array float>: tableau des coordonnées avec les colonnes ['X', 'Y', 'Xt_amont', 'Xt_aval']
        """
        new_coord = np.empty(len(coord), dtype=self.points.dtype)
        # FIXME: avoid copying in using np.lib.recfunctions.append_fields?
        for var in ['X', 'Y', 'xt', 'Xt_amont', 'Xt_aval']:  # copy existing columns
            new_coord[var] = coord[var]
        new_coord['zone'] = index_zone
        new_coord['lit'] = index_lit
        new_coord['Xl'] = self.profils_travers[index_zone].dist_proj_axe * (1 - xl) + \
                          self.profils_travers[index_zone + 1].dist_proj_axe * xl
        new_coord['xl'] = xl
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

    def export_triangulation_dict(self):
        """
        @brief: Export triangulation with vertices in cartesian coordinates
        """
        return {'vertices': np.array(np.column_stack((self.points['X'], self.points['Y']))),
                'segments': self.segments}

    def export_floworiented_triangulation_dict(self):
        """
        @brief: Export triangulation with vertices in flow-oriented coordinates
        """
        u = self.points['Xt_amont'] * (1 - self.points['xl']) + self.points['Xt_aval'] * self.points['xl']
        v = self.points['Xl']
        return {'vertices': np.array(np.column_stack((u, v))),
                'segments': self.segments}

    def build_initial_profiles(self):
        logger.info("~> Interpolation sur les profils existants en prenant en compte "
                    "le passage des lignes de contraintes")
        for i in range(len(self.profils_travers)):
            cur_profil = self.profils_travers[i]
            logger.debug(cur_profil)

            # Recherche des limites communes "amont"
            if i == 0:
                common_limites_id_1 = cur_profil.limites.keys()
            else:
                common_limites_id_1 = cur_profil.common_limits(self.profils_travers[i - 1].limites.keys())

            # Recherche des limites communes "aval"
            if i == len(self.profils_travers) - 1:
                common_limites_id_2 = cur_profil.limites.keys()
            else:
                common_limites_id_2 = cur_profil.common_limits(self.profils_travers[i + 1].limites.keys())

            # Union ordonnée des limites amont/aval
            limites_id = cur_profil.common_limits(list(set(common_limites_id_1).union(common_limites_id_2)))

            first_lit = True
            for j, (id1, id2) in enumerate(zip(limites_id, limites_id[1:])):
                lit = cur_profil.extraire_lit(id1, id2)
                coord_int = lit.interp_coord_along_lit_auto(self.pas_trans, self.nb_pts_trans)

                if first_lit:
                    cur_profil.get_limit_by_id(id1)['id_pt'] = self.i_pt + 1
                else:
                    coord_int = coord_int[1:]

                if i == 0:
                    coord_int = rename_fields(coord_int, {'Xt': 'Xt_amont'})
                    coord_int = append_fields(coord_int, 'Xt_aval', np.zeros(len(coord_int)), usemask=False)
                    self.add_points(coord_int, 0.0, i, j)
                else:
                    coord_int = rename_fields(coord_int, {'Xt': 'Xt_aval'})
                    coord_int = append_fields(coord_int, 'Xt_amont', np.zeros(len(coord_int)), usemask=False)
                    self.add_points(coord_int, 1.0, i - 1, j)

                cur_profil.get_limit_by_id(id2)['id_pt'] = self.i_pt

                # Ajoute les nouveaux segments
                new_i_pt = np.arange(cur_profil.get_limit_by_id(id1)['id_pt'],
                                     cur_profil.get_limit_by_id(id2)['id_pt'] + 1)
                self.add_segments_from_node_list(new_i_pt)

                if first_lit:
                    first_lit = False

    def build_interp(self, lignes_contraintes, pas_long, constant_ech_long):
        """
        Build interpolation, add points and segments

        @param lignes_contraintes <[LigneContrainte]>: lignes de contrainte
        @param pas_long <float>
        @param constant_ech_long <bool>
        """
        self.build_initial_profiles()

        ### BOUCLE SUR L'ESPACE INTER-PROFIL
        logger.info("~> Construction du maillage par zone interprofils puis par lit")

        for i, (prev_profil, next_profil) in enumerate(zip(self.profils_travers, self.profils_travers[1:])):
            logger.debug("> Zone n°{} : entre {} et {}".format(i, prev_profil, next_profil))

            if constant_ech_long:
                nb_pts_inter = prev_profil.calcul_nb_pts_inter(next_profil, pas_long)
                Xp_adm_list = np.linspace(0.0, 1.0, num=nb_pts_inter + 2)[1:-1]

            # Recherche des limites communes entre les deux profils
            common_limites_id = prev_profil.common_limits(next_profil.limites.keys())
            logger.debug("Limites de lits communes : {}".format(list(common_limites_id)))

            if len(common_limites_id) < 2:
                raise TatooineException("Aucune interpolation pour l'intervalle {}, entre {} et {} ({} limites communes)".format(
                    i, prev_profil, next_profil, len(common_limites_id)))

            else:
                first_lit = True
                ### BOUCLE SUR LES LITS (= MORCEAU(X) DE PROFIL)
                for j, (id1, id2) in enumerate(zip(common_limites_id, common_limites_id[1:])):
                    pt_list_L1 = []
                    pt_list_L2 = []

                    logger.debug("Lit {}-{}".format(id1, id2))

                    # Extraction d'une partie des profils
                    lit_1 = prev_profil.extraire_lit(id1, id2)
                    lit_2 = next_profil.extraire_lit(id1, id2)

                    # Abscisses curvilignes le long de la ligne de contrainte
                    (Xp_profil1_L1, Xp_profil1_L2) = prev_profil.get_Xt_lignes(id1, id2)
                    (Xp_profil2_L1, Xp_profil2_L2) = next_profil.get_Xt_lignes(id1, id2)
                    dXp_L1 = Xp_profil2_L1 - Xp_profil1_L1
                    dXp_L2 = Xp_profil2_L2 - Xp_profil1_L2

                    if dXp_L1 < 0:
                        raise TatooineException("La ligne {} n'est pas orientée dans le même ordre que les profils"
                                                .format(id1))
                    if dXp_L2 < 0:
                        raise TatooineException("La ligne {} n'est pas orientée dans le même ordre que les profils"
                                                .format(id2))

                    if not constant_ech_long:
                        nb_pts_inter = math.ceil(min(dXp_L1, dXp_L2)/pas_long) - 1
                        Xp_adm_list = np.linspace(0.0, 1.0, num=nb_pts_inter + 2)[1:-1]

                    L1_coord_int = lignes_contraintes[id1].coord_sampling_along_line(Xp_profil1_L1, Xp_profil2_L1,
                                                                                     Xp_adm_list)
                    L2_coord_int = lignes_contraintes[id2].coord_sampling_along_line(Xp_profil1_L2, Xp_profil2_L2,
                                                                                     Xp_adm_list)

                    ### BOUCLE SUR LES LIGNES
                    for k in range(nb_pts_inter):
                        Xp = Xp_adm_list[k]
                        P1 = Point(tuple(L1_coord_int[k]))
                        P2 = Point(tuple(L2_coord_int[k]))

                        if self.nb_pts_trans is None:
                            nb_pts_trans = math.ceil(P1.distance(P2) / self.pas_trans) + 1
                        else:
                            nb_pts_trans = self.nb_pts_trans
                        array = lit_1.interp_coord_linear(lit_2, Xp, nb_pts_trans)
                        lit_int = Lit(array, ['Xt', 'xt'])
                        lit_int.move_between_targets(P1, P2)
                        coord_int = lit_int.array[['X', 'Y', 'xt', 'Xt_amont', 'Xt_aval']]  # Ignore `Xt`
                        pt_list_L1.append(self.i_pt + 1)

                        if not first_lit:
                            # ignore le 1er point car la ligne de contrainte a déjà été traitée
                            coord_int = coord_int[1:]

                        self.add_points(coord_int, Xp, i, j)

                        pt_list_L2.append(self.i_pt)

                    pt_list_L2 = np.array([prev_profil.get_limit_by_id(id2)['id_pt']] + pt_list_L2 +
                                          [next_profil.get_limit_by_id(id2)['id_pt']])
                    self.add_segments_from_node_list(pt_list_L2)

                    if first_lit:
                        pt_list_L1 = np.array([prev_profil.get_limit_by_id(id1)['id_pt']] + pt_list_L1 +
                                              [next_profil.get_limit_by_id(id1)['id_pt']])
                        self.add_segments_from_node_list(pt_list_L1)
                        first_lit = False

    def corr_bathy_on_epis(self, epis, dist_corr_epi):
        raise NotImplementedError
        logger.info("~> Correction de la bathymétrie autour des épis")
        if self.var_names() != ['Z']:
            raise TatooineException("Impossible de corriger les épis sur les couches sédimentaires")
        for epi in epis:
            epi_geom = epi.coord.convert_as_linestring()
            for i, coord in enumerate(self.points):
                pt_node = Point((coord['X'], coord['Y']))
                if epi_geom.distance(pt_node) < dist_corr_epi:
                    Xt_proj = epi_geom.project(pt_node)
                    pt_proj = epi_geom.interpolate(Xt_proj)
                    epi.coord.values['Z'][i] = pt_proj.z

    def build_mesh(self, in_floworiented_crs=False):
        logger.info("~> Calcul du maillage")
        if in_floworiented_crs:
            tri = self.export_floworiented_triangulation_dict()
            self.triangle = triangle.triangulate(tri, opts='p')
            self.triangle['vertices'] = self.export_triangulation_dict()['vertices']  # overwrite by cartesian coordinates
        else:
            tri = self.export_triangulation_dict()
            self.triangle = triangle.triangulate(tri, opts='p')
        if len(self.points) != len(self.triangle['vertices']):
            if len(self.points) < len(self.triangle['vertices']):
                logger.error("New nodes are:")
                ori_points = np.column_stack((self.points['X'], self.points['Y']))
                ori_combined = ori_points[:, 0] * ori_points[:, 1] / (ori_points[:, 0] + ori_points[:, 1])
                new_points = self.triangle['vertices']
                new_combined = new_points[:, 0] * new_points[:, 1] / (new_points[:, 0] + new_points[:, 1])
                diff = np.setxor1d(ori_combined, new_combined)
                logger.error(new_points[np.isin(new_combined, diff)])
            raise TatooineException("Mesh is corrupted... %i vs %i nodes." % (
                len(self.points), len(self.triangle['vertices'])))
        if 'triangles' not in self.triangle:
            raise TatooineException("Mesh was not generated, no triangle found!")
        logger.info(self.summary())

    def summary(self):
        try:
            nnode, nelem = len(self.triangle['vertices']), len(self.triangle['triangles'])
        except KeyError:
            raise TatooineException("La génération du maillage a échouée!")
        return "Génération d'un maillage avec {} noeuds et {} éléments".format(nnode, nelem)

    def export_points(self, path):
        if path.endswith('.xyz'):
            logger.info("~> Exports en xyz des points")
            with open(path, 'wb') as fileout:
                z_array = self.interp_values_from_profiles()[0, :]
                np.savetxt(fileout, np.vstack((self.points['X'], self.points['Y'], z_array)).T,
                           delimiter=' ', fmt='%.{}f'.format(DIGITS))
        elif path.endswith('.shp'):
            logger.info("~> Exports en shp des points")
            z_array = self.interp_values_from_profiles()[0, :]
            with shapefile.Writer(path, shapeType=shapefile.POINT) as w:
                w.field('zone', 'N', decimal=6)
                w.field('lit', 'N', decimal=6)
                w.field('Xt_amont', 'N', decimal=6)
                w.field('Xt_aval', 'N', decimal=6)
                w.field('xt', 'N', decimal=6)
                w.field('xl', 'N', decimal=6)
                w.field('Z', 'N', decimal=6)
                for row, z in zip(self.points, z_array):
                    w.point(row['X'], row['Y'])
                    w.record(**{'zone': float(row['zone']), 'lit': float(row['lit']),
                                'Xt_amont': row['Xt_amont'], 'Xt_aval': row['Xt_aval'], 'xt': row['xt'],
                                'xl': row['xl'], 'Z': z})
        else:
            raise NotImplementedError("Seuls les formats shp et xyz sont supportés pour les semis de points")

    def export_segments(self, path):
        if path.endswith('.shp'):
            logger.info("~> Exports en shp des segments")
            with shapefile.Writer(path, shapeType=shapefile.POLYLINE) as w:
                w.field('id_seg', 'N', decimal=6)
                for i, (node1, node2) in enumerate(self.segments):
                    point1 = self.points[node1]
                    point2 = self.points[node2]
                    w.line([[[point1['X'], point1['Y']], [point2['X'], point2['Y']]]])
                    w.record(id_seg=i)
        else:
            raise NotImplementedError("Seul le format shp est supporté pour les segments")

    def export_profiles(self, path):
        """
        /!\ Pas cohérent si constant_ech_long est différent de True
        """
        values = self.interp_values_from_profiles()
        if path.endswith('.georefC'):
            with open(path, 'w') as out_geo:
                for dist in np.unique(self.points['Xl']):
                    pos = self.points['Xl'] == dist
                    points = self.points[pos]

                    # Compute Xt  (FIXME: rather keep from previous calculations...)
                    Xt = np.sqrt(np.power(np.ediff1d(points['X'], to_begin=0.), 2) +
                                 np.power(np.ediff1d(points['Y'], to_begin=0.), 2))
                    Xt = Xt.cumsum()
                    points = append_fields(points, 'Xt', Xt, usemask=False)

                    for i, row in enumerate(points):
                        if i == 0:
                            positions_str = ' %f %f %f %f' % (row['X'], row['Y'], points[-1]['X'], points[-1]['Y'])
                            positions_str += ' AXE %f %f' % (row['X'], row['Y'])  # FIXME: not the axis position...
                            out_geo.write('Profil Bief_0 %s %f%s\n' % ('P' + str(dist), dist, positions_str))

                        layers_str = ' ' + ' '.join([COURLIS_FLOAT_FMT % x for x in values[:, pos][:, i]])

                        out_geo.write('%f%s B %f %f\n' % (row['Xt'], layers_str, row['X'], row['Y']))
            return

        lines = []
        for dist in np.unique(self.points['Xl']):
            pos = self.points['Xl'] == dist
            line = geometry.Polyline([(x, y, z) for (x, y), z in zip(self.points[pos][['X', 'Y']], values[0, :])])
            line.add_attribute(dist)
            lines.append(line)

        if path.endswith('.i3s'):
            with bk.Write(path) as out_i3s:
                out_i3s.write_header()
                out_i3s.write_lines(lines, [l.attributes()[0] for l in lines])
        elif path.endswith('.shp'):
            shp.write_shp_lines(path, shapefile.POLYLINEZ, lines, 'Z')
        else:
            raise NotImplementedError("Seuls les formats i3s, georefC et shp (POLYLINEZ) sont supportés pour écrire le "
                                      "fichier de profils en travers")

    def export_mesh(self, path):
        """TODO: export multiple variables in t3s and LandXML"""
        logger.info("~> Écriture du maillage")

        nnode, nelem = len(self.triangle['vertices']), len(self.triangle['triangles'])
        if path.endswith('.t3s'):
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
:WrittenBy                tatooinemesher
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
                np.savetxt(fileout, np.column_stack((self.triangle['vertices'], self.interp_values_from_profiles()[0, :])),
                           delimiter=' ', fmt='%.{}f'.format(DIGITS))

                # Tableau des éléments (connectivité)
                np.savetxt(fileout, self.triangle['triangles'] + 1, delimiter=' ', fmt='%i')

        elif path.endswith('.xml'):
            env = Environment(
                loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'data')))
            template = env.get_template("LandXML_template.xml")
            template_render = template.render(
                nodes=np.round(np.column_stack((self.triangle['vertices'], self.interp_values_from_profiles()[0, :])),
                               DIGITS),
                ikle=self.triangle['triangles'] + 1
            )

            # Écriture du fichier XML
            with open(path, 'w') as fileout:
                fileout.write(template_render)

        elif path.endswith('.slf'):

            with Serafin.Write(path, LANG, overwrite=True) as resout:
                output_header = Serafin.SerafinHeader(title='%s (written by tatooinemesher)' % os.path.basename(path),
                                                      lang=LANG)
                output_header.from_triangulation(self.triangle['vertices'], self.triangle['triangles'] + 1)

                for var in self.var_names():
                    if var in basic_2D_vars_IDs:
                        output_header.add_variable_from_ID(var)
                    else:
                        output_header.add_variable_str(var, var, '')
                resout.write_header(output_header)

                resout.write_entire_frame(output_header, 0.0, self.interp_values_from_profiles())

        else:
            raise NotImplementedError("Seuls les formats t3d, xml et slf sont supportés pour les maillages")

    def interp_1d_values_from_profiles(self):
        """Interpolate values in 1D (lateral + longitudinal) from profiles"""
        new_values = np.zeros((self.profils_travers[0].coord.nb_var(), len(self.points)))
        for i_zone in np.unique(self.points['zone']):
            filter_points = self.points['zone'] == i_zone
            section_us = self.profils_travers[i_zone]
            section_ds = self.profils_travers[i_zone + 1]
            xt_us = section_us.coord.array['Xt']
            xt_ds = section_ds.coord.array['Xt']
            xt_us_target = self.points['Xt_amont'][filter_points]
            xt_ds_target = self.points['Xt_aval'][filter_points]

            for i, var in enumerate(self.var_names()):
                values_us = section_us.coord.values[var]
                values_ds = section_ds.coord.values[var]

                if self.interp_values == 'LINEAR':
                    new_values_us = np.interp(xt_us_target, xt_us, values_us)
                    new_values_ds = np.interp(xt_ds_target, xt_ds, values_ds)

                elif self.interp_values == 'B-SPLINE':
                    splrep_us = interpolate.splrep(xt_us, values_us)
                    splrep_ds = interpolate.splrep(xt_ds, values_ds)
                    new_values_us = interpolate.splev(xt_us_target, splrep_us)
                    new_values_ds = interpolate.splev(xt_ds_target, splrep_ds)

                elif self.interp_values == 'AKIMA':
                    new_values_us = interpolate.Akima1DInterpolator(xt_us, values_us)(xt_us_target)
                    new_values_ds = interpolate.Akima1DInterpolator(xt_ds, values_ds)(xt_ds_target)

                elif self.interp_values == 'PCHIP':
                    new_values_us = interpolate.pchip_interpolate(xt_us, values_us, xt_us_target)
                    new_values_ds = interpolate.pchip_interpolate(xt_ds, values_ds, xt_ds_target)

                elif self.interp_values == 'CUBIC_SPLINE':
                    new_values_us = interpolate.CubicSpline(xt_us, values_us)(xt_us_target)
                    new_values_ds = interpolate.CubicSpline(xt_ds, values_ds)(xt_ds_target)

                else:
                    raise NotImplementedError

                new_values[i, filter_points] = new_values_us * (1 - self.points['xl'][filter_points]) + \
                                               new_values_ds * self.points['xl'][filter_points]
        return new_values

    def interp_2d_values_from_profiles(self):
        """Interpolate values in 2D from profiles"""
        ux = np.array([], dtype=np.float)
        vy = np.array([], dtype=np.float)
        new_xt = self.points['xt']
        new_xl = self.points['xl'] + self.points['zone']
        new_values = np.zeros((self.profils_travers[0].coord.nb_var(), len(self.points)))
        for i, profile in enumerate(self.profils_travers):
            first_xt = profile.get_limit_by_idx(0)['Xt_profil']
            last_xt = profile.get_limit_by_idx(-1)['Xt_profil']
            xt = (profile.coord.array['Xt'] - first_xt)/(last_xt - first_xt)
            ux = np.concatenate((ux, xt))
            vy = np.concatenate((vy, np.array([i] * profile.nb_points)))

        for j, var in enumerate(self.var_names()):
            z = np.array([], dtype=np.float)
            for profile in self.profils_travers:
                z = np.concatenate((z, profile.coord.values[var]))

            if self.interp_values == 'BIVARIATE_SPLINE':
                interp_bivariate_spline = interpolate.SmoothBivariateSpline(ux, vy, z, kx=3, ky=3)
                for k, (u, v) in enumerate(zip(new_xt, new_xl)):
                    new_values[j, k] = interp_bivariate_spline(u, v)[0][0]

            else:
                if self.interp_values == 'BILINEAR':
                    method = 'linear'
                elif self.interp_values == 'BICUBIC':
                    method = 'cubic'
                else:
                    raise NotImplementedError
                new_values[j, :] = interpolate.griddata((ux, vy), z, (new_xt, new_xl), method=method)

        return new_values

    def interp_values_from_profiles(self):
        if self.interp_values in ('BILINEAR', 'BICUBIC', 'BIVARIATE_SPLINE'):
            return self.interp_2d_values_from_profiles()
        else:
            return self.interp_1d_values_from_profiles()

    def interp_from_values_at_profiles(self, values_at_profiles):
        """Interpolate values from profiles"""
        nb_var = values_at_profiles.shape[1]
        values = np.zeros((nb_var, len(self.points)))
        xl_all = self.points['zone'] + self.points['xl']
        for i_var in range(nb_var):
            values[i_var, :] = np.interp(xl_all, np.arange(len(self.profils_travers), dtype=np.float),
                                         values_at_profiles[:, i_var])
        return values

    @staticmethod
    def get_merge_triangulation(mesh_constr_list):
        """
        Merge multiple distinct triangulations (no merge at boundaries or overlaps).
        @param mesh_constr_list <[MeshConstructor]>: list of mesh constructors
        """
        out_tri = {}
        for mesh_constr in mesh_constr_list:
            if not out_tri:  # set initial values from first iteration
                out_tri['triangles'] = mesh_constr.triangle['triangles']
                out_tri['vertices'] = mesh_constr.triangle['vertices']
            else:  # concatenate with current sub-mesh for next iterations
                out_tri['triangles'] = np.vstack((out_tri['triangles'],
                                                  out_tri['vertices'].shape[0] + mesh_constr.triangle['triangles']))
                out_tri['vertices'] = np.vstack((out_tri['vertices'], mesh_constr.triangle['vertices']))
        return out_tri


class GeometryRequestException(Exception):
    """Custom exception for geometry parser"""
    def __init__(self, message):
        super().__init__(message)
