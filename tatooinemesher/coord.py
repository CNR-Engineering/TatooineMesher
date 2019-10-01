from copy import deepcopy
import numpy as np
from numpy.lib.recfunctions import append_fields
from shapely.geometry import LineString, Point

from tatooinemesher.utils import logger, strictly_increasing, TatooineException


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
