from collections import OrderedDict
from copy import deepcopy
import math
import matplotlib.pyplot as plt
import numpy as np
import os.path
from pyteltools.geom import BlueKenue as bk, Shapefile as shp
import shapefile
from shapely.geometry import LineString, MultiPoint, Point

from tatooinemesher.coord import Coord
from tatooinemesher.utils import float_vars, get_field_index, get_intersections, logger, \
    strictly_increasing, TatooineException


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
