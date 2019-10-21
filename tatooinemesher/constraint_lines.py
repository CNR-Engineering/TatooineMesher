import numpy as np
from pyteltools.geom import BlueKenue as bk, Shapefile as shp
import shapefile
from shapely.geometry import LineString

from tatooinemesher.interp.cubic_hermite_spline import CubicHermiteSpline
from tatooinemesher.utils import float_vars, TatooineException


class ConstraintLine:
    """
    ConstraintLine: open polyline used to separate beds

    ### Attributes
    - id <integer>:  unique identifier (automatic numbering starting from 0)
    - nb_points <int>: number of points
    - coord <2D-array float>: coordinates (X and Y)
    - geom <shapely.geometry.LineString>: 2D geometry
    - interp: coordinates interpolator

    ### Methods
    - get_lines_from_file
    - get_lines_and_set_limits_from_sections
    - build_interp_linear
    - build_interp_chs
    - coord_sampling_along_line
    """
    def __init__(self, id, coord, interp_coord='LINEAR'):
        """
        Create a ConstraintLine instance from X and Y coordinates
        @id <integer>: unique identifier (automatic numbering starting from 0)
        @coord <Coord>: coordinates (X and Y)
            Ex: [(x1, y1), (x2, y2), ..., (xn, yn)]
        @interp_coord <str>: coordinates interpolation method (among: 'LINEAR', 'CARDINAL' and 'FINITE_DIFF')
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
        return "ConstraintLine #{} ({} points)".format(self.id, self.nb_points)

    @staticmethod
    def get_lines_from_file(filename, interp_coord='LINEAR'):
        """
        Returns a list of ConstraintLine from an input file
        TODO 1: Value is ignored in i2s file format
        """
        lines = []
        if filename is not None:
            if filename.endswith('.i2s'):
                with bk.Read(filename) as in_i2s:
                    in_i2s.read_header()
                    for i, line in enumerate(in_i2s.get_open_polylines()):
                        lines.append(ConstraintLine(i, list(line.polyline().coords), interp_coord))

            elif filename.endswith('.shp'):
                if shp.get_shape_type(filename) not in (shapefile.POLYLINE, shapefile.POLYLINEZ, shapefile.POLYLINEM):
                    raise TatooineException("The type of file %s is not POLYLINEZ[M]" % filename)
                for i, line in enumerate(shp.get_open_polylines(filename)):
                    lines.append(ConstraintLine(i, list(line.polyline().coords), interp_coord))

            else:
                raise NotImplementedError("Only shp and i3s formats are supported for constraint lines")

        return lines

    @staticmethod
    def get_lines_and_set_limits_from_sections(section_seq, interp_coord='LINEAR'):
        """
        @brief: Returns a list of ConstraintLine from an sequence of cross-sections
        @param section_seq <CrossSectionSequence>: sequence of cross-sections
        @param interp_coord <str>: interpolation method
        """
        # Build 2 constraint lines from cross-section bounds
        first_coords = []
        last_coords = []
        for section in section_seq:
            first_coords.append(section.geom.coords[0][:2])
            last_coords.append(section.geom.coords[-1][:2])
        lines = [ConstraintLine(0, first_coords, interp_coord),
                 ConstraintLine(1, last_coords, interp_coord)]

        # Set limits
        for line_id, line in enumerate(lines):
            for section, Xt_line in zip(section_seq, line.Xt):
                Xt_section = section.coord.array['Xt'][0] if line_id == 0 else section.coord.array['Xt'][-1]
                intersection = section.geom.interpolate(Xt_section)
                section.add_limit(line_id, Xt_section, Xt_line, intersection)

        return lines

    def build_interp_linear(self):
        """
        @brief: Build a double linear interpolator for X and Y coordinates
        """
        def interp_xy_linear(Xt_new):
            coord_int = []
            for dist in Xt_new:
                point = self.geom.interpolate(dist)
                coord_int.append(point.coords[0][:2])
            return np.array(coord_int, dtype=float_vars(['X', 'Y']))
        return interp_xy_linear

    def build_interp_chs(self, tan_method):
        """
        @brief: Build a Cubic Hermite Spline interpolator for X and Y coordinates
        """
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
        @brief: Extract coordinates interpolated along line between Xp1 and Xp2
            at requested dimensionless curvilinear distances
        @param: Xp_adm_int <1D-array>: dimensionless curvilinear distances between Xp1 and Xp2 (0=Xp1 and 1=Xp2)
        @param: Xp1 <float>: starting curvilinear distance
        @param: Xp2 <float>: ending curvilinear distance
        """
        # Building list of curvilinear distance in meters
        Xp = (1 - Xp_adm_int)*Xp1 + Xp_adm_int*Xp2
        # Use coordinate interpolator
        return self.interp(Xp)
