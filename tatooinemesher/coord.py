from copy import deepcopy
import numpy as np
from numpy.lib.recfunctions import append_fields
from shapely.geometry import LineString, Point

from tatooinemesher.utils import logger, strictly_increasing, TatooineException


class Coord:
    """
    Coord: consecutive 2D or 3D points (e.g. polylines)

    ### Attributes
    - array: structured array with coordinates ('X', 'Y'), variables (see coord_labels)
        and eventually distance(s) ('Xt', 'xt')
    - coord_labels: list of variable names (TOCHECK)
    - values: structured array with values (different variables are possible)

    ### Methods
    - nb_var
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
                logger.warn("Duplicated points are removed")
                # Suppression des doublons (points superposés dans la polyligne)
                points2keep = np.ediff1d(self.array['Xt'], to_begin=1.) != 0.
                self.array = self.array[points2keep]

    def nb_var(self):
        if self.values is None:
            return 0
        else:
            return len(self.values.dtype.names)

    def compute_Xt(self):
        """
        Compute cumulated curvilinear distance `Xt`
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
        Compute dimensionless curvilinear distance `xt` (from 0 to 1)
        /!\ Column `Xt` has to exist
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
        @brief: Shift current instance to fit bounds to p1 and p2
            (by a linear weighted translation from p1 to p2)
        @param p1 <Point>: starting point
        @param p2 <Point>: ending point
        """
        array = deepcopy(self.array)
        xt = array['xt']
        self.array['X'] = array['X'] + (1-xt)*(p1.x-array['X'][0]) + xt*(p2.x-array['X'][-1])
        self.array['Y'] = array['Y'] + (1-xt)*(p1.y-array['Y'][0]) + xt*(p2.y-array['Y'][-1])

    def compute_xp(self):
        """
        Compute dimensionless from starting to ending point distance projetée adimensionnée sur droite début->fin
        """
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
