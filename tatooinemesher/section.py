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


class CrossSection:
    """
    CrossSection: represents a cross-section with possibly multiple variables at points

    ### Attributes
    - id <integer|str>: unique identifier
    - label <str>: type of section (`Cross-section` or `Transverse constraint line`)
    - coord <Coord>: coordinates (X and Y) with all variables
    - nb_points <int>: number of points
    - geom <shapely.geometry.LineString>: 2D geometry
    - limits <OrderedDict>: ordered dictionary with `id_ligne` as keys and a tuple with
        (Xt_section, Xt_line, X, Y, ...) as values
    - dist_proj_axe (used for <CrossSectionSequence>)

    ### Methods
    - add_limit
    - get_limit_by_id
    - get_limit_by_idx
    - get_Xt_lines
    - find_and_add_limit
    - common_limits
    - extract_bed
    - sort_limits
    - compute_nb_pts_inter
    - project_straight_line
    - get_segments
    - get_angles
    - export_plot_crosswise
    """
    def __init__(self, id, coord, label='Cross-section'):
        """
        Create a cross-section from coordinates X and Y
        Limits are not set at the object instantiation

        @param id <integer|str>: unique identifier
        @label <str>: type of section (`Cross-section` or `Transverse constraint line`)
        @param coords <[tuple]>: sequence of X and Y coordinates
        """
        self.id = id
        self.label = label
        self.coord = Coord(np.array(coord, dtype=float_vars(['X', 'Y'])), ['Xt', 'xt'])
        self.nb_points = len(self.coord.array)
        self.geom = LineString(coord)  # FIXME: might contains duplicated points

        self.limits = OrderedDict()
        self.dist_proj_axe = -1

    def __repr__(self):
        return "{} #{} ({} points)".format(self.label, self.id, self.nb_points)

    def add_limit(self, line_id, Xt_section, Xt_line, point):
        """Add a new limit"""
        z_values = {}
        for label in self.coord.coord_labels:
            z_values[label] = np.interp(Xt_section, self.coord.array['Xt'], self.coord.array[label])
        self.limits[line_id] = {'Xt_section': Xt_section, 'Xt_line': Xt_line,
                                'X': point.x, 'Y': point.y, **z_values}

    def get_limit_by_id(self, line_id):
        return self.limits[line_id]

    def get_limit_by_idx(self, idx):
        return self.limits[list(self.limits.keys())[idx]]

    def get_Xt_lines(self, id1, id2):
        return self.get_limit_by_id(id1)['Xt_line'], self.get_limit_by_id(id2)['Xt_line']

    def find_and_add_limit(self, constraint_line, dist_max=None):
        """
        @param constraint_line <shapely.geometry.LineString>: 2D constraint line
        @param dist_max <float>: maximum search distance to rescue intersections for limits
        """
        if self.geom.intersects(constraint_line.geom):
            intersection = self.geom.intersection(constraint_line.geom)

            if isinstance(intersection, MultiPoint):
                logger.warn("Intersection between '{}' and '{}' contains multiple points, "
                            "only the first is kept.".format(self, constraint_line))
                intersection = intersection[0]

            if isinstance(intersection, Point):
                # Compute projections
                Xt_section = self.geom.project(intersection)
                Xt_line = constraint_line.geom.project(intersection)
                self.add_limit(constraint_line.id, Xt_section, Xt_line, intersection)

            else:
                raise TatooineException("Intersection between '{}' and '{}' is empty or not supported: {}".format(
                    self, constraint_line, type(intersection)))
        else:
            if dist_max is not None:
                distance = self.geom.distance(constraint_line.geom)
                if distance < dist_max:
                    for i, coord in enumerate(constraint_line.coord):
                        # Try to find a point of the constraint line which is in the vicinity of the current section
                        point = Point(coord)
                        dist = self.geom.distance(point)
                        if dist < dist_max:
                            # A point is found and is considered
                            Xt_line = constraint_line.geom.project(point)
                            Xt_section = self.geom.project(point)
                            intersection = self.geom.interpolate(Xt_section)
                            self.add_limit(constraint_line.id, Xt_section, Xt_line, intersection)
                            logger.debug("Add a limit with the line {} after {} iterations (distance = {})"
                                         .format(constraint_line.id, i, dist))
                            break

    def common_limits(self, limit_ids):
        """
        @brief: List common limits between 2 cross-sections
        @param limit_ids: list of limit identifiers of other cross-section
        """
        out_limits = []
        for limit_id in self.limits.keys():
            if limit_id in limit_ids:
                out_limits.append(limit_id)
        return out_limits

    def extract_bed(self, bed1_id, bed2_id):
        """
        @brief: Extract coordinates of a portion of a cross-section between 2 limits
            /!\ bed1_id and bed2_id should be "ordered" correctly, otherwise an exception is raised
        @return <Bed>: structured array with columns ('X', 'Y', 'Xt', 'xt')
        """
        limit1 = self.get_limit_by_id(bed1_id)
        limit2 = self.get_limit_by_id(bed2_id)

        Xt1 = limit1['Xt_section']
        Xt2 = limit2['Xt_section']

        # Check that Xt are increasing from bed1_id to bed2_id
        if Xt1 > Xt2:
            raise TatooineException("Order of beds {} and {} leads to decreasing Xt values for {}".format(
                bed1_id, bed2_id, self))

        Xt_section = self.coord.array['Xt']
        sub_coord = self.coord.array[np.logical_and(Xt_section >= Xt1, Xt_section <= Xt2)]

        # Add starting point if necessary
        if Xt1 not in Xt_section:
            row = np.array([tuple(limit1[var] if var not in ('Xt', 'xt') else Xt1 for var in sub_coord.dtype.names)],
                           dtype=sub_coord.dtype)
            row['xt'] = 0.0
            sub_coord = np.insert(sub_coord, 0, row)

        # Add last point if necessary
        if Xt2 not in Xt_section:
            row = np.array([tuple(limit2[var] if var not in ('Xt', 'xt') else Xt2 for var in sub_coord.dtype.names)],
                           dtype=sub_coord.dtype)
            row['xt'] = 1.0
            sub_coord = np.append(sub_coord, row)

        # Check order of points
        if not strictly_increasing(sub_coord['Xt']):
            logger.debug("/!\ Xt values are not strictly increasing")  # FIXME: It should not happen!
            logger.debug(sub_coord['Xt'])
            logger.debug("Please check the following limits below:")
            logger.debug(limit1)
            logger.debug(limit2)
            points_to_keep = np.ediff1d(sub_coord['Xt'], to_begin=1.) != 0.
            sub_coord = sub_coord[points_to_keep]

        return Bed(sub_coord, ['xt'])

    def sort_limits(self):
        """
        @brief: Sort limits by increasing Xt values
        """
        self.limits = OrderedDict(sorted(self.limits.items(), key=lambda x: x[1]['Xt_section']))

    def compute_nb_pts_inter(self, other, long_dist):
        """
        @brief: Compute the number of intermediate cross-sections required
        /!\ This counter includes points at bounds (start and end)
        @param other <CrossSection>: other cross-section
        @param long_dist <float>: longitudinal space step (in m)
        """
        dist_min_profil = self.geom.distance(other.geom)
        return math.ceil(dist_min_profil / long_dist) - 1

    def project_straight_line(self):
        """
        @brief: Planar projection on a straight line joining first and last point of the cross-section
        """
        if self.limits:
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

    def export_plot_crosswise(self, fig_path, overwrite=False):
        x = self.coord.array['Xt']
        z = self.coord.values['Z']

        fig, ax1 = plt.subplots(figsize=(16, 9))

        ax1.set_xlabel('Xt (m)')

        color = 'tab:red'
        ax1.set_ylabel('Z', color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.plot(x, z, marker='o', color=color, label='Z')

        color = 'tab:green'
        for limit_name, limit in self.limits.items():
            ax1.axvline(x=limit['Xt_section'], linestyle='-', color=color)

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


class Bed(Coord):
    """
    Bed: portion of a cross-section

    ### Attributes
    - coord

    ### Methods
    - interp_coord_along_bed_auto
    - interp_coord_along_bed
    - interp_coord_linear
    """
    def interp_coord_along_bed_auto(self, lat_step, nb_pts_lat=None):
        if nb_pts_lat is None:
            nb_pts_lat = math.ceil((self.array['Xt'][-1] - self.array['Xt'][0])/lat_step) + 1
        Xt_adm_list = np.linspace(0., 1., num=nb_pts_lat)
        return self.interp_coord_along_bed(Xt_adm_list)

    def interp_coord_along_bed(self, Xt_adm_list):
        """
        @brief: Interpolate coordinates along bed at requested dimensionless curvilinear distances
        @param Xt_adm_list <1D-array float>: dimensionless distances (ranging from 0 and 1)
        """
        nb_pts_lat = len(Xt_adm_list)
        array = np.empty(nb_pts_lat, dtype=float_vars(Coord.XY + ['Xt', 'xt']))
        for label in Coord.XY + ['Xt']:
            array[label] = np.interp(Xt_adm_list, self.array['xt'], self.array[label])
        array['xt'] = Xt_adm_list
        return array

    def interp_coord_linear(self, other, coeff, nb_pts_lat):
        """
        @brief: Linear interpolation between 2 beds
        @param other <Bed>: portion of the other cross-section (upstream or downstream)
        @param coeff <float>: weight between self and other (0=self, 1=other)
        @param nb_pts_lat <int>: number of points crosswise
        """
        Xt_adm_list = np.linspace(0., 1., num=nb_pts_lat)

        array_1 = self.interp_coord_along_bed(Xt_adm_list)
        array_2 = other.interp_coord_along_bed(Xt_adm_list)

        array = np.empty(nb_pts_lat, dtype=float_vars(Coord.XY + ['xt', 'Xt_upstream', 'Xt_downstream']))
        for var in Coord.XY:
            array[var] = (1 - coeff)*array_1[var] + coeff*array_2[var]
        array['Xt_upstream'] = array_1['Xt']
        array['Xt_downstream'] = array_2['Xt']
        return array


class CrossSectionSequence:
    """
    CrossSectionSequence: ordered list of cross-sections

    ### Attributes
    - section_list <[CrossSection]>: ordered list of cross-sections

    ### Methods
    - add_section
    - from_file
    - find_and_add_limits
    - check_intersections
    - compute_dist_proj_axe
    - sort_by_dist
    - export_sections_shp
    - add_constant_layer
    """

    def __init__(self):
        self.section_list = []

    def add_section(self, section):
        self.section_list.append(section)

    @staticmethod
    def from_file(filename, label, field_id=None, project_straight_line=False):
        section_seq = CrossSectionSequence()

        if filename.endswith('.i3s'):
            with bk.Read(filename) as in_i3s:
                in_i3s.read_header()
                for i, line in enumerate(in_i3s.get_open_polylines()):
                    line_id = i if field_id is None else line.attributes()[0]  # Use `Value` if field is not None
                    z_array = np.array([(coord[2],) for coord in line.polyline().coords], dtype=float_vars('Z'))
                    line = line.to_2d()
                    section = CrossSection(line_id, list(line.polyline().coords), label)
                    section.coord.values = z_array
                    section_seq.add_section(section)

        elif filename.endswith('.shp'):
            shp_type = shp.get_shape_type(filename)
            if shp_type in (shapefile.POLYLINEZ, shapefile.POLYLINEM):
                field_id_index = get_field_index(filename, field_id)
                for i, line in enumerate(shp.get_open_polylines(filename)):
                    line_id = i if field_id is None else line.attributes()[field_id_index]
                    z_array = np.array([(coord[2],) for coord in line.polyline().coords], dtype=float_vars(['Z']))
                    line = line.to_2d()
                    section = CrossSection(line_id, list(line.polyline().coords), label)
                    section.coord.values = z_array
                    section_seq.add_section(section)

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
                        section = CrossSection(last_point_id, coords, label)
                        section.coord.values = z_array
                        section_seq.add_section(section)
                        coords, z_layers = [], []
                    coords.append(point[:2])
                    z_layers.append((point[2],) + tuple(attributes[index] for index in field_indexes))
                    last_point_id = point_id
                z_array = np.array(z_layers, dtype=float_vars(['Z'] + field_names))
                section = CrossSection(last_point_id, coords, label)
                section.coord.values = z_array
                section_seq.add_section(section)

            else:
                raise TatooineException("The type of file %s is not POINTZ or POLYLINEZ[M]" % filename)

        else:
            raise NotImplementedError("Only shp and i3s formats are supported for cross-sections")

        if project_straight_line:
            for section in section_seq:
                section.project_straight_line()
        return section_seq

    def __add__(self, other):
        new_suite = deepcopy(self)
        new_suite.section_list = self.section_list + other.suite
        return new_suite

    def __getitem__(self, index):
        if isinstance(index, slice):
            new_suite = deepcopy(self)
            new_suite.section_list = self.section_list[index]
            return new_suite
        else:
            return self.section_list[index]

    def __len__(self):
        return len(self.section_list)

    def __repr__(self):
        return [section for section in self.section_list].__repr__()

    def find_and_add_limits(self, constraint_lines, dist_max):
        """
        @param constraint_lines <[ConstraintLine]>: list of constraint lines
        @param dist_max <float>: maximum search distance to rescue intersections for limits
        """
        logger.info("~> Looking for limits")
        for i, section in enumerate(self):
            for constraint_line in constraint_lines:
                section.find_and_add_limit(constraint_line, dist_max)
            section.sort_limits()

            limits = section.limits.keys()  # only to print
            logger.debug("> {}".format(section))
            logger.debug("{} limits found with lines {}".format(len(limits), list(limits)))

    def check_intersections(self):
        """Display warning with intersections details"""
        logger.info("~> Checking that cross-sections do not intersect")
        intersections = get_intersections([section.geom for section in self])
        if intersections:
            logger.warn("Following limits are found:")
            for (i, j) in intersections:
                logger.warn("- between '{}' and '{}'".format(self[i], self[j]))

    def compute_dist_proj_axe(self, axe_geom, dist_max):
        """
        @brief: Compute distance along hydraulic axis
        @param axe_geom <shapely.geometry.LineString>: hydraulic axis (/!\ Beware of its orientation)
        @param dist_max <float>: maximum search distance to rescue intersections for limits
        """
        logger.info("~> Compute distances along hydraulic axis to order cross-sections")
        to_keep_list = []
        for section in self:
            section_geom = section.geom
            if section_geom.intersects(axe_geom):
                intersection = section_geom.intersection(axe_geom)
                if isinstance(intersection, Point):
                    section.dist_proj_axe = axe_geom.project(intersection)
                else:
                    raise TatooineException("Intersection between '{}' and the hydraulic axis "
                                            "is not a unique point".format(section))
            else:
                if dist_max is not None:
                    for pos in (0, -1):
                        dist = section_geom.distance(Point(axe_geom.coords[pos]))
                        if dist < dist_max:
                            section.dist_proj_axe = 0.0 if pos == 0 else axe_geom.length
            if section.dist_proj_axe == -1:
                logger.warn("{} do not intersect the hydraulic axis (distance = {}m) and is ignored".format(
                    section, section.geom.distance(axe_geom)))
                to_keep_list.append(False)
            else:
                to_keep_list.append(True)

        self.section_list = [p for p, to_keep in zip(self.section_list, to_keep_list) if to_keep]

    def sort_by_dist(self):
        self.section_list = sorted(self.section_list, key=lambda x: x.dist_proj_axe)

    def export_sections_shp(self, out_path):
        """
        Write a shapefile with section as 3D LineString
        @param out_path <str>: output file name
        """
        with shapefile.Writer(out_path, shapeType=shapefile.POLYLINEZ) as w:
            w.field('profil_id', 'C')
            for profil in self.section_list:
                array = profil.coord.array
                z_array = profil.coord.values['Z']
                coords = [(row['X'], row['Y'], z) for row, z in zip(array, z_array)]
                w.linez([coords])
                w.record(profil_id=str(profil.id))

    def add_constant_layer(self, name, thickness):
        for section in self.section_list:
            coord = section.coord
            coord.add_single_layer(name, coord.array['Z'] - thickness)
