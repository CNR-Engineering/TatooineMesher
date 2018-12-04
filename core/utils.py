import logging
import numpy as np
import shapefile
import sys

from pyteltools.geom import BlueKenue as bk, Shapefile as shp


logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(handler)
logger.setLevel(logging.INFO)


def get_intersections(linestrings):
    """
    @brief: Chercher les intersections entre toutes les lignes
    @param linestrings <[shapely.geometry.LineString]>: lignes à analyser
    """
    intersections = []
    for i, l1 in enumerate(linestrings):
        for j, l2 in enumerate(linestrings[i+1:]):
            if l1.intersects(l2):
                intersections.append((i, j+i+1))
    return intersections


def float_vars(varname_list):
    return [(label, np.float) for label in varname_list]


def strictly_increasing(array):
    """
    @brief: Check if array is sorted
    @param array <1D-array float>: array to check
    """
    return all(x < y for x, y in zip(array, array[1:]))


def get_axe_hydraulique(infile_axe):
    """
    @brief: Extraction de l'unique ligne du fichier i2s d'entrée
    @param infile_axe <str>: chemin vers le fichier contenant l'axe hydraulique
    @return <shapely.geometry.LineString>: polyligne représentant l'axe
    """
    if infile_axe.endswith('.i2s'):
        with bk.Read(infile_axe) as in_i2s:
            in_i2s.read_header()
            lines = list(in_i2s.get_open_polylines())
    elif infile_axe.endswith('.shp'):
        if shp.get_shape_type(infile_axe) not in (shapefile.POLYLINE, shapefile.POLYLINEZ, shapefile.POLYLINEM):
            raise TatooineException("Le fichier %s n'est pas de type POLYLINEZ" % infile_axe)
        lines = list(shp.get_open_polylines(infile_axe))
    else:
        raise NotImplementedError("Seuls les formats i2s et shp sont supportés pour l'axe hydraulique")
    nb_lines = len(lines)
    if nb_lines != 1:
        raise TatooineException("Le fichier '{}' contient {} polylignes au lieu d'une seule pour définir "
                                "l'axe hydraulique".format(infile_axe, nb_lines))
    return lines[0].polyline()


class TatooineException(Exception):
    """!
    @brief Custom exception for TatooineMesher
    """
    def __init__(self, message):
        """!
        @param message <str>: error message description
        """
        super().__init__(message)
        self.message = message
