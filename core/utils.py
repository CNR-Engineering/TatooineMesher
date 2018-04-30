import numpy as np
import sys

from pyteltools.geom import BlueKenue as bk


def get_intersections(linestrings):
    """
    @brief: Chercher les intersections entre toutes les lignes
    @param linestrings <[LineString]>: lignes à analyser
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
    """Extraction de l'unique ligne du fichier i2s d'entrée"""
    with bk.Read(infile_axe) as in_i2s:
        in_i2s.read_header()
        lines = [line.polyline() for line in in_i2s.get_open_polylines()]
    nb_lines = len(lines)
    if nb_lines != 1:
        sys.exit("Le fichier '{}' contient {} polylignes au lieu d'une seule pour définir l'axe hydraulique".format(infile_axe, nb_lines))
    return lines[0]