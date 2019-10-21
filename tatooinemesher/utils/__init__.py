import logging
from math import ceil, sqrt
import numpy as np
import shapefile

from crue10.utils import logger as crue10_logger
from pyteltools.geom import BlueKenue as bk, Shapefile as shp
from pyteltools.utils.log import set_logger_level as set_pyteltools_logger_level


logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(handler)


def get_intersections(linestrings):
    """
    @brief: Find intersections between all lines
    @param linestrings <[shapely.geometry.LineString]>: list of lines to analyze
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


def get_hydraulic_axis(infile_axis):
    """
    @brief: Extract a unique line from i2s input file
    @param infile_axis <str>: path to file
    @return <shapely.geometry.LineString>: polyline representing the hydraulic axis
    """
    if infile_axis.endswith('.i2s'):
        with bk.Read(infile_axis) as in_i2s:
            in_i2s.read_header()
            lines = list(in_i2s.get_open_polylines())
    elif infile_axis.endswith('.shp'):
        if shp.get_shape_type(infile_axis) not in (shapefile.POLYLINE, shapefile.POLYLINEZ, shapefile.POLYLINEM):
            raise TatooineException("The type of file %s is not POLYLINE[ZM]" % infile_axis)
        lines = list(shp.get_open_polylines(infile_axis))
    else:
        raise NotImplementedError("Only shp and i2s formats are supported for hydraulic axis")
    nb_lines = len(lines)
    if nb_lines != 1:
        raise TatooineException("The file '{}' contains {} polylines instead of a unique line to define "
                                "the hydraulic axis".format(infile_axis, nb_lines))
    return lines[0].polyline()


def resample_2d_line(coord, dist_max):
    """
    @brief: resamples a 2D polyline by preserving its initial points
    @param coord <[tuple]>: vertices coordinates as list of tuples [(x1, y1), (x2, y2), ...]
    @param dist_max <float>: maximal distance to refine segments
    """
    new_coord = [coord[0]]  # add first point

    for coord_line in zip(coord, coord[1:]):
        # Iterate over each segment
        [(x1, y1), (x2, y2)] = coord_line
        length = sqrt((x1 - x2)**2 + (y1 - y2)**2)

        if length > dist_max:
            # Nombre de points pour la nouvelle ligne (en incluant les 2 points d'origine)
            nb_pts = ceil(length/dist_max) + 1
            # Sampling with prescribe number of points
            #   (and ignore first point which was in the previous
            xnew = np.linspace(x1, x2, num=nb_pts)[1:]
            ynew = np.linspace(y1, y2, num=nb_pts)[1:]

            for x, y in zip(xnew, ynew):
                new_coord.append((x, y))

        else:
            # Add ending point
            new_coord.append((x2, y2))

    return new_coord


def get_field_index(filename, field_id):
    if field_id is not None:
        names, _ = shp.get_attribute_names(filename)
        try:
            return names.index(field_id)
        except ValueError:
            raise TatooineException("The field `%s` does not exist" % field_id)


def set_logger_level(set_to_debug):
    level = logging.DEBUG if set_to_debug else logging.INFO
    logger.setLevel(level)
    crue10_logger.setLevel(level)
    set_pyteltools_logger_level(level)


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
