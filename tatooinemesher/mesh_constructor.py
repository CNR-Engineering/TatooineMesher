from jinja2 import Environment, FileSystemLoader
import math
import numpy as np
from numpy.lib.recfunctions import append_fields, rename_fields
import os.path
from pyteltools.geom import BlueKenue as bk, Shapefile as shp
from pyteltools.geom import geometry
from pyteltools.slf import Serafin
from pyteltools.slf.variable.variables_2d import basic_2D_vars_IDs
from scipy import interpolate
import shapefile
from shapely.geometry import Point
import time
import triangle

from tatooinemesher.section import Bed
from tatooinemesher.utils import float_vars, logger, TatooineException


DIGITS = 4  # for csv and xml exports
COURLIS_FLOAT_FMT = "%.6f"


class MeshConstructor:
    """
    MeshConstructor: construire un maillage et interpoler les valeurs

    ### Attributs
    - section_seq <CrossSectionSequence>
    - lat_step <float>: pas transversal (m)
    - nb_pts_lat <int>: nombre de noeuds transversalement
    - interp_values <str>: interpolation method
    - nb_var <int>: number of variables
    - points: structured array with columns ['X', 'Y', 'Xt_upstream', 'Xt_downstream', 'Xl', 'xl', 'zone', 'bed']
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
    - export_sections
    - export_mesh
    - interp_values_from_geom
    - get_merge_triangulation
    """

    POINTS_DTYPE = float_vars(["X", "Y", "xt", "Xt_upstream", "Xt_downstream", "Xl", "xl"]) + [
        (var, int) for var in ("zone", "bed")
    ]
    POINTS_FP_DTYPE = float_vars(["X", "Y", "Z"])

    def __init__(self, section_seq=[], lat_step=None, nb_pts_lat=None, interp_values="LINEAR"):
        self.section_seq = section_seq
        self.lat_step = lat_step
        self.nb_pts_lat = nb_pts_lat
        self.interp_values = interp_values
        if self.section_seq:
            self.nb_var = self.section_seq[0].coord.nb_var()
        else:
            self.nb_var = 0

        self.points = np.empty(0, dtype=MeshConstructor.POINTS_DTYPE)
        self.nodes_values = np.empty([0, 0, 0], dtype=float)
        self.i_pt = int(-1)
        self.segments = np.empty([0, 2], dtype=int)
        self.triangle = {}  # filled by `build_mesh`

        self.casiers_nodes_idx = []
        self.nodes_fp = np.empty(0, dtype=MeshConstructor.POINTS_FP_DTYPE)

    @property
    def nb_nodes_in_riverbed(self):
        return len(self.points)

    @property
    def nb_nodes(self):
        return self.nb_nodes_in_riverbed + len(self.nodes_fp)

    @property
    def has_floodplain(self):
        return len(self.casiers_nodes_idx) != 0

    def var_names(self):
        return list(self.section_seq[0].coord.values.dtype.names)

    def add_points(self, coord, zone_index, xl, bed_index):
        """!
        @brief: Add vertices/nodes
        @param coord <2D-array float>: table with columns 'X', 'Y', 'Xt_upstream' and 'Xt_downstream'
        """
        if self.casiers_nodes_idx:
            raise TatooineException("Impossible to add points in river bed after having considered the floodplain")

        new_coord = np.empty(len(coord), dtype=self.points.dtype)
        # FIXME: avoid copying in using np.lib.recfunctions.append_fields?
        for var in ["X", "Y", "xt", "Xt_upstream", "Xt_downstream"]:  # copy existing columns
            new_coord[var] = coord[var]
        new_coord["zone"] = zone_index
        new_coord["bed"] = bed_index
        new_coord["Xl"] = (
            self.section_seq[zone_index].dist_proj_axe * (1 - xl) + self.section_seq[zone_index + 1].dist_proj_axe * xl
        )
        new_coord["xl"] = xl
        self.i_pt += len(new_coord)
        self.points = np.hstack((self.points, new_coord))

    def add_floodplain_mesh(self, triangulation, points):
        pos_start = self.nb_nodes
        pos_end = pos_start + len(points)
        self.nodes_fp = np.hstack((self.nodes_fp, points))
        self.casiers_nodes_idx.append((pos_start, pos_end))
        self._append_triangulation(triangulation)

    def add_segments(self, segments):
        """!
        @brief: Add segments to force from an array of node indexes
        @param segments <2D-array int>: array with node indexes
        """
        self.segments = np.vstack((self.segments, segments))

    def add_segments_from_node_list(self, node_list):
        """
        @brief: Add segments to force from an list of consecutive node indexes
        @param node_list <1D-array int>: list of node indexes
        """
        new_segments = np.column_stack((node_list[:-1], node_list[1:]))
        self.add_segments(new_segments)

    def _append_triangulation(self, triangulation):
        """
        @brief: Add a submesh from a triangulation
        @param triangulation <dict>: dictionary with `triangles` and `vertices` keys
        """
        if not self.triangle:
            self.triangle["triangles"] = triangulation["triangles"]
            self.triangle["vertices"] = triangulation["vertices"]
        else:
            nb_elem = self.triangle["vertices"].shape[0]
            self.triangle["triangles"] = np.vstack((self.triangle["triangles"], nb_elem + triangulation["triangles"]))
            self.triangle["vertices"] = np.vstack((self.triangle["vertices"], triangulation["vertices"]))

    def append_mesh_constr(self, mesh_constr):
        """
        @brief: Append a local mesh to the current instance (adjacent or superimposed meshes are not merged)
        @param mesh_constr <MeshConstructor>: MeshConstructor to combine with current instance
        """
        if self.nb_nodes == 0:
            self.nb_var = mesh_constr.nb_var
        else:
            if self.nb_var != mesh_constr.nb_var:
                raise RuntimeError

        self.section_seq += mesh_constr.section_seq

        points = mesh_constr.points
        if len(self.points) != 0:
            last_zone = self.points["zone"].max()
            points["zone"] += last_zone + 2
        self.points = np.hstack((self.points, points))

        self._append_triangulation(mesh_constr.triangle)

    def export_triangulation_dict(self):
        """
        @brief: Export triangulation with nodes in geometric coordinates
        """
        return {"vertices": np.array(np.column_stack((self.points["X"], self.points["Y"]))), "segments": self.segments}

    def export_floworiented_triangulation_dict(self):
        """
        @brief: Export triangulation with nodes in flow-oriented coordinates
        """
        u = self.points["Xt_upstream"] * (1 - self.points["xl"]) + self.points["Xt_downstream"] * self.points["xl"]
        v = self.points["Xl"]
        return {"vertices": np.array(np.column_stack((u, v))), "segments": self.segments}

    def build_initial_profiles(self):
        logger.info("~> Interpolate initial cross-sections taking into account limits")
        for i in range(len(self.section_seq)):
            curr_profile = self.section_seq[i]
            logger.debug(curr_profile)

            # Looking for upstream common limits
            if i == 0:
                common_limits_id_1 = curr_profile.limits.keys()
            else:
                common_limits_id_1 = curr_profile.common_limits(self.section_seq[i - 1].limits.keys())

            # Looking for downstream common limits
            if i == len(self.section_seq) - 1:
                common_limits_id_2 = curr_profile.limits.keys()
            else:
                common_limits_id_2 = curr_profile.common_limits(self.section_seq[i + 1].limits.keys())

            # Ordered union of upstream and downstream limits
            limits_id = curr_profile.common_limits(list(set(common_limits_id_1).union(common_limits_id_2)))

            first_bed = True
            for j, (id1, id2) in enumerate(zip(limits_id, limits_id[1:])):
                bed = curr_profile.extract_bed(id1, id2)
                coord_int = bed.interp_coord_along_bed_auto(self.lat_step, self.nb_pts_lat)

                if first_bed:
                    curr_profile.get_limit_by_id(id1)["id_pt"] = self.i_pt + 1
                else:
                    coord_int = coord_int[1:]

                if i == 0:
                    coord_int = rename_fields(coord_int, {"Xt": "Xt_upstream"})
                    coord_int = append_fields(coord_int, "Xt_downstream", np.zeros(len(coord_int)), usemask=False)
                    self.add_points(coord_int, i, 0.0, j)
                else:
                    coord_int = rename_fields(coord_int, {"Xt": "Xt_downstream"})
                    coord_int = append_fields(coord_int, "Xt_upstream", np.zeros(len(coord_int)), usemask=False)
                    self.add_points(coord_int, i - 1, 1.0, j)

                curr_profile.get_limit_by_id(id2)["id_pt"] = self.i_pt

                # Add new segments
                new_i_pt = np.arange(
                    curr_profile.get_limit_by_id(id1)["id_pt"], curr_profile.get_limit_by_id(id2)["id_pt"] + 1
                )
                self.add_segments_from_node_list(new_i_pt)

                if first_bed:
                    first_bed = False

    def build_interp(self, constraint_lines, long_step, constant_long_disc):
        """
        Build interpolation, add points and segments

        @param constraint_lines <[ConstraintLine]>: list of constraint lines
        @param long_step <float>: longitudinal space step
        @param constant_long_disc <bool>
        """
        nb_pts_inter = 0
        self.build_initial_profiles()

        # LOOP ON ZONES (between 2 consecutive cross-sections)
        logger.info("~> Building mesh per zone and then per bed")

        for i, (prev_section, next_section) in enumerate(zip(self.section_seq, self.section_seq[1:])):
            logger.debug("> Zone n°{} : between {} and {}".format(i, prev_section, next_section))

            if constant_long_disc:
                nb_pts_inter = prev_section.compute_nb_pts_inter(next_section, long_step)
                Xp_adm_list = np.linspace(0.0, 1.0, num=nb_pts_inter + 2)[1:-1]

            # Looking for common limits between cross-sections
            common_limits_id = prev_section.common_limits(next_section.limits.keys())
            logger.debug("Common limits: {}".format(list(common_limits_id)))

            if len(common_limits_id) < 2:
                raise TatooineException(
                    "No interpolation in the interval %i, between %s and %s (%i common limits)"
                    % (i, prev_section, next_section, len(common_limits_id))
                )

            else:
                first_bed = True
                # LOOP ON BEDS
                for j, (id1, id2) in enumerate(zip(common_limits_id, common_limits_id[1:])):
                    pt_list_L1 = []
                    pt_list_L2 = []

                    logger.debug("Bed {}-{}".format(id1, id2))

                    # Extraction of cross-section portions (= beds)
                    bed_1 = prev_section.extract_bed(id1, id2)
                    bed_2 = next_section.extract_bed(id1, id2)

                    # Curvilinear abscissa along constraint lines
                    (Xp_profil1_L1, Xp_profil1_L2) = prev_section.get_Xt_lines(id1, id2)
                    (Xp_profil2_L1, Xp_profil2_L2) = next_section.get_Xt_lines(id1, id2)
                    dXp_L1 = Xp_profil2_L1 - Xp_profil1_L1
                    dXp_L2 = Xp_profil2_L2 - Xp_profil1_L2

                    if dXp_L1 < 0:
                        raise TatooineException("The constraint line {} is not oriented correctly".format(id1))
                    if dXp_L2 < 0:
                        raise TatooineException("The constraint line {} is not oriented correctly".format(id2))

                    if not constant_long_disc:
                        nb_pts_inter = math.ceil(min(dXp_L1, dXp_L2) / long_step) - 1
                        Xp_adm_list = np.linspace(0.0, 1.0, num=nb_pts_inter + 2)[1:-1]

                    L1_coord_int = constraint_lines[id1].coord_sampling_along_line(
                        Xp_profil1_L1, Xp_profil2_L1, Xp_adm_list
                    )
                    L2_coord_int = constraint_lines[id2].coord_sampling_along_line(
                        Xp_profil1_L2, Xp_profil2_L2, Xp_adm_list
                    )

                    # LOOP ON INTERMEDIATE CROSS-SECTIONS
                    for k in range(nb_pts_inter):
                        Xp = Xp_adm_list[k]
                        P1 = Point(tuple(L1_coord_int[k]))
                        P2 = Point(tuple(L2_coord_int[k]))

                        if self.nb_pts_lat is None:
                            nb_pts_lat = math.ceil(P1.distance(P2) / self.lat_step) + 1
                        else:
                            nb_pts_lat = self.nb_pts_lat
                        array = bed_1.interp_coord_linear(bed_2, Xp, nb_pts_lat)
                        bed_int = Bed(array, ["Xt", "xt"])
                        bed_int.move_between_targets(P1, P2)
                        coord_int = bed_int.array[["X", "Y", "xt", "Xt_upstream", "Xt_downstream"]]  # Ignore `Xt`
                        pt_list_L1.append(self.i_pt + 1)

                        if not first_bed:
                            # ignore first point because the constraint line was already considered
                            coord_int = coord_int[1:]

                        self.add_points(coord_int, i, Xp, j)

                        pt_list_L2.append(self.i_pt)

                    pt_list_L2 = np.array(
                        [prev_section.get_limit_by_id(id2)["id_pt"]]
                        + pt_list_L2
                        + [next_section.get_limit_by_id(id2)["id_pt"]]
                    )
                    self.add_segments_from_node_list(pt_list_L2)

                    if first_bed:
                        pt_list_L1 = np.array(
                            [prev_section.get_limit_by_id(id1)["id_pt"]]
                            + pt_list_L1
                            + [next_section.get_limit_by_id(id1)["id_pt"]]
                        )
                        self.add_segments_from_node_list(pt_list_L1)
                        first_bed = False

    def corr_bathy_on_epis(self, epis, dist_corr_epi):
        raise NotImplementedError  # TODO
        logger.info("~> Correction bathymetry around 'épis'")
        if self.var_names() != ["Z"]:
            raise TatooineException("Impossible to correct elevation of 'épis' for the sediment layers")
        for epi in epis:
            epi_geom = epi.coord.convert_as_linestring()
            for i, coord in enumerate(self.points):
                pt_node = Point((coord["X"], coord["Y"]))
                if epi_geom.distance(pt_node) < dist_corr_epi:
                    Xt_proj = epi_geom.project(pt_node)
                    pt_proj = epi_geom.interpolate(Xt_proj)
                    epi.coord.values["Z"][i] = pt_proj.z

    def build_mesh(self, in_floworiented_crs=False, opts="p"):
        """
        @brief: Build mesh under constraints
        @param in_floworiented_crs <bool>: boolean to which coordinate system is usec
        @param opts <str>: options for the triangulation.
            `p` - Triangulates a Planar Straight Line Graph.
            `q` - Quality mesh generation with no angles smaller than 20 degrees.
                  An alternate minimum angle may be specified after the `q`.
            `a` - Imposes a maximum triangle area constraint. A fixed area constraint (that applies to every triangle)
                  may be specified after the `a`, or varying areas may be read from the input dictionary.
        """
        logger.info("~> Building mesh")
        if in_floworiented_crs:
            tri = self.export_floworiented_triangulation_dict()
            self.triangle = triangle.triangulate(tri, opts=opts)
            self.triangle["vertices"] = self.export_triangulation_dict()[
                "vertices"
            ]  # overwrite by cartesian coordinates
        else:
            tri = self.export_triangulation_dict()
            self.triangle = triangle.triangulate(tri, opts=opts)

        if opts == "p":  # Check that vertices correspond to points
            if len(self.points) != len(self.triangle["vertices"]):
                if len(self.points) < len(self.triangle["vertices"]):
                    logger.error("New nodes are:")
                    ori_points = np.column_stack((self.points["X"], self.points["Y"]))
                    ori_combined = ori_points[:, 0] * ori_points[:, 1] / (ori_points[:, 0] + ori_points[:, 1])
                    new_points = self.triangle["vertices"]
                    new_combined = new_points[:, 0] * new_points[:, 1] / (new_points[:, 0] + new_points[:, 1])
                    diff = np.setxor1d(ori_combined, new_combined)
                    logger.error(new_points[np.isin(new_combined, diff)])
                raise TatooineException(
                    "Mesh is corrupted... %i vs %i nodes." % (len(self.points), len(self.triangle["vertices"]))
                )

        if "triangles" not in self.triangle:
            raise TatooineException("Mesh was not generated, no triangle found!")
        logger.info(self.summary())

    def summary(self):
        try:
            nnode, nelem = len(self.triangle["vertices"]), len(self.triangle["triangles"])
        except KeyError:
            raise TatooineException("The generation of the mesh failed!")
        return "Mesh with {} nodes and {} elements".format(nnode, nelem)

    def export_points(self, path):
        if path.endswith(".xyz"):
            logger.info("~> Exports points in xyz")
            with open(path, "wb") as fileout:
                z_array = self.interp_values_from_geom()[0, :]
                np.savetxt(
                    fileout,
                    np.vstack((self.points["X"], self.points["Y"], z_array)).T,
                    delimiter=" ",
                    fmt="%.{}f".format(DIGITS),
                )

        elif path.endswith(".shp"):
            logger.info("~> Exports points in shp")
            z_array = self.interp_values_from_geom()[0, :]
            with shapefile.Writer(path, shapeType=shapefile.POINT) as w:
                w.field("zone", "N", decimal=6)
                w.field("bed", "N", decimal=6)
                w.field("Xt_upstream", "N", decimal=6)
                w.field("Xt_downstream", "N", decimal=6)
                w.field("xt", "N", decimal=6)
                w.field("xl", "N", decimal=6)
                w.field("Z", "N", decimal=6)
                for row, z in zip(self.points, z_array):
                    w.point(row["X"], row["Y"])
                    w.record(
                        **{
                            "zone": float(row["zone"]),
                            "bed": float(row["bed"]),
                            "Xt_upstream": row["Xt_upstream"],
                            "Xt_downstream": row["Xt_downstream"],
                            "xt": row["xt"],
                            "xl": row["xl"],
                            "Z": z,
                        }
                    )

        else:
            raise TatooineException("Only shp and xyz are supported for points set")

    def export_segments(self, path):
        if path.endswith(".shp"):
            logger.info("~> Export segments in a shp file")
            with shapefile.Writer(path, shapeType=shapefile.POLYLINE) as w:
                w.field("id_seg", "N", decimal=6)
                for i, (node1, node2) in enumerate(self.segments):
                    point1 = self.points[node1]
                    point2 = self.points[node2]
                    w.line([[[point1["X"], point1["Y"]], [point2["X"], point2["Y"]]]])
                    w.record(id_seg=i)
        else:
            raise TatooineException("Only the shp format is supported for segments")

    def export_sections(self, path):
        """
        Export generated profiles in a shp, i3s or georefC file
        /!\ Not relevant if constant_long_disc is False
        TODO: Use class MascaretGeoFile
        """
        values = self.interp_values_from_geom()
        if path.endswith(".georefC"):
            with open(path, "w") as out_geo:
                for dist in np.unique(self.points["Xl"]):
                    pos = self.points["Xl"] == dist
                    points = self.points[pos]

                    # Compute Xt  (FIXME: rather keep from previous calculations...)
                    Xt = np.sqrt(
                        np.power(np.ediff1d(points["X"], to_begin=0.0), 2)
                        + np.power(np.ediff1d(points["Y"], to_begin=0.0), 2)
                    )
                    Xt = Xt.cumsum()
                    points = append_fields(points, "Xt", Xt, usemask=False)

                    for i, row in enumerate(points):
                        if i == 0:
                            positions_str = " %f %f %f %f" % (row["X"], row["Y"], points[-1]["X"], points[-1]["Y"])
                            positions_str += " AXE %f %f" % (row["X"], row["Y"])  # FIXME: not the axis position...
                            out_geo.write("Profil Bief_0 %s %f%s\n" % ("P" + str(dist), dist, positions_str))

                        layers_str = " " + " ".join([COURLIS_FLOAT_FMT % x for x in values[:, pos][:, i]])

                        out_geo.write("%f%s B %f %f\n" % (row["Xt"], layers_str, row["X"], row["Y"]))
            return

        lines = []
        for dist in np.unique(self.points["Xl"]):
            pos = self.points["Xl"] == dist
            line = geometry.Polyline([(x, y, z) for (x, y), z in zip(self.points[pos][["X", "Y"]], values[0, :])])
            line.add_attribute(dist)
            lines.append(line)

        if path.endswith(".i3s"):
            with bk.Write(path) as out_i3s:
                out_i3s.write_header()
                out_i3s.write_lines(lines, [l.attributes()[0] for l in lines])

        elif path.endswith(".shp"):
            shp.write_shp_lines(path, shapefile.POLYLINEZ, lines, "Z")

        else:
            raise NotImplementedError(
                "Only the shp (POLYLINEZ), i3s and georefC formats are supported for "
                "the generated cross-sections file"
            )

    def export_mesh(self, path, lang="en"):
        """
        Export generated mesh in slf, t3s or LandXML
        TODO: export multiple variables in t3s and LandXML
        """
        logger.info("~> Write generated mesh")

        nnode, nelem = len(self.triangle["vertices"]), len(self.triangle["triangles"])
        if path.endswith(".t3s"):
            with open(path, "w", newline="") as fileout:
                # Write header
                date = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
                fileout.write(
                    """#########################################################################
:FileType t3s  ASCII  EnSim 1.0
# Canadian Hydraulics Centre/National Research Council (c) 1998-2012
# DataType                 2D T3 Scalar Mesh
#
:Application              BlueKenue
:Version                  3.3.4
:WrittenBy                TatooineMesher
:CreationDate             {}
#
#------------------------------------------------------------------------
#
:NodeCount {}
:ElementCount {}
:ElementType  T3
#
:EndHeader
""".format(
                        date, nnode, nelem
                    )
                )

            with open(path, mode="ab") as fileout:
                # Table with x, y, z coordinates
                np.savetxt(
                    fileout,
                    np.column_stack((self.triangle["vertices"], self.interp_values_from_geom()[0, :])),
                    delimiter=" ",
                    fmt="%.{}f".format(DIGITS),
                )

                # Table with elements (connectivity)
                np.savetxt(fileout, self.triangle["triangles"] + 1, delimiter=" ", fmt="%i")

        elif path.endswith(".xml"):
            env = Environment(
                loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.realpath(__file__)), "data"))
            )
            template = env.get_template("LandXML_template.xml")
            template_render = template.render(
                nodes=np.round(
                    np.column_stack((self.triangle["vertices"], self.interp_values_from_geom()[0, :])), DIGITS
                ),
                ikle=self.triangle["triangles"] + 1,
            )

            # Write XML file
            with open(path, "w") as fileout:
                fileout.write(template_render)

        elif path.endswith(".slf"):

            with Serafin.Write(path, lang, overwrite=True) as resout:
                output_header = Serafin.SerafinHeader(
                    title="%s (Written by TatooineMesher)" % os.path.basename(path), lang=lang
                )
                output_header.from_triangulation(self.triangle["vertices"], self.triangle["triangles"] + 1)

                for var in self.var_names():
                    if var in basic_2D_vars_IDs:
                        output_header.add_variable_from_ID(var)
                    else:
                        output_header.add_variable_str(var, var, "")
                resout.write_header(output_header)

                resout.write_entire_frame(output_header, 0.0, self.interp_values_from_geom())

        else:
            raise NotImplementedError("Only slf, t3s and xml formats are supported for the output mesh")

    def compute_values_in_floodplain(self):
        """Fill bottom elevation and set nan for other variables"""
        values = np.empty((self.nb_var, len(self.nodes_fp)))
        values.fill(np.nan)
        values[0, :] = self.nodes_fp["Z"]
        return values

    def interp_1d_values_from_profiles(self):
        """Interpolate values in 1D (lateral + longitudinal) from profiles"""
        new_values = np.zeros((self.nb_var, self.nb_nodes_in_riverbed))
        for i_zone in np.unique(self.points["zone"]):
            filter_points = self.points["zone"] == i_zone
            section_us = self.section_seq[i_zone]
            section_ds = self.section_seq[i_zone + 1]
            xt_us = section_us.coord.array["Xt"]
            xt_ds = section_ds.coord.array["Xt"]
            xt_us_target = self.points["Xt_upstream"][filter_points]
            xt_ds_target = self.points["Xt_downstream"][filter_points]

            for i, var in enumerate(self.var_names()):
                values_us = section_us.coord.values[var]
                values_ds = section_ds.coord.values[var]

                if self.interp_values == "LINEAR":
                    new_values_us = np.interp(xt_us_target, xt_us, values_us)
                    new_values_ds = np.interp(xt_ds_target, xt_ds, values_ds)

                elif self.interp_values == "B-SPLINE":
                    splrep_us = interpolate.splrep(xt_us, values_us)
                    splrep_ds = interpolate.splrep(xt_ds, values_ds)
                    new_values_us = interpolate.splev(xt_us_target, splrep_us)
                    new_values_ds = interpolate.splev(xt_ds_target, splrep_ds)

                elif self.interp_values == "AKIMA":
                    new_values_us = interpolate.Akima1DInterpolator(xt_us, values_us)(xt_us_target)
                    new_values_ds = interpolate.Akima1DInterpolator(xt_ds, values_ds)(xt_ds_target)

                elif self.interp_values == "PCHIP":
                    new_values_us = interpolate.pchip_interpolate(xt_us, values_us, xt_us_target)
                    new_values_ds = interpolate.pchip_interpolate(xt_ds, values_ds, xt_ds_target)

                elif self.interp_values == "CUBIC_SPLINE":
                    new_values_us = interpolate.CubicSpline(xt_us, values_us)(xt_us_target)
                    new_values_ds = interpolate.CubicSpline(xt_ds, values_ds)(xt_ds_target)

                else:
                    raise NotImplementedError

                new_values[i, filter_points] = (
                    new_values_us * (1 - self.points["xl"][filter_points])
                    + new_values_ds * self.points["xl"][filter_points]
                )
        return new_values

    def interp_2d_values_from_profiles(self):
        """Interpolate values in 2D from profiles"""
        ux = np.array([], dtype=float)
        vy = np.array([], dtype=float)
        new_xt = self.points["xt"]
        new_xl = self.points["xl"] + self.points["zone"]
        new_values = np.zeros((self.nb_var, len(self.points)))
        for i, profile in enumerate(self.section_seq):
            first_xt = profile.get_limit_by_idx(0)["Xt_profil"]
            last_xt = profile.get_limit_by_idx(-1)["Xt_profil"]
            xt = (profile.coord.array["Xt"] - first_xt) / (last_xt - first_xt)
            ux = np.concatenate((ux, xt))
            vy = np.concatenate((vy, np.array([i] * profile.nb_points)))

        for j, var in enumerate(self.var_names()):
            z = np.array([], dtype=float)
            for profile in self.section_seq:
                z = np.concatenate((z, profile.coord.values[var]))

            if self.interp_values == "BIVARIATE_SPLINE":
                interp_bivariate_spline = interpolate.SmoothBivariateSpline(ux, vy, z, kx=3, ky=3)
                for k, (u, v) in enumerate(zip(new_xt, new_xl)):
                    new_values[j, k] = interp_bivariate_spline(u, v)[0][0]

            else:
                if self.interp_values == "BILINEAR":
                    method = "linear"
                elif self.interp_values == "BICUBIC":
                    method = "cubic"
                else:
                    raise NotImplementedError
                new_values[j, :] = interpolate.griddata((ux, vy), z, (new_xt, new_xl), method=method)

        return new_values

    def interp_values_from_geom(self):
        """Interpolate values from geometry at cross-sections and casiers"""
        values = np.empty((self.nb_var, self.nb_nodes))
        values.fill(np.nan)
        if self.interp_values in ("BILINEAR", "BICUBIC", "BIVARIATE_SPLINE"):
            values[:, : self.nb_nodes_in_riverbed] = self.interp_2d_values_from_profiles()
        else:
            values[:, : self.nb_nodes_in_riverbed] = self.interp_1d_values_from_profiles()
        values[:, self.nb_nodes_in_riverbed : self.nb_nodes] = self.compute_values_in_floodplain()
        return values

    def interp_values_from_res(self, values_at_sections, z_at_casiers, pos_z):
        """Interpolate values from results at cross-sections and casiers"""
        nb_var = values_at_sections.shape[1]
        values = np.empty((nb_var, self.nb_nodes))
        values.fill(np.nan)
        xl_all = self.points["zone"] + self.points["xl"]
        for i_var in range(nb_var):
            values[i_var, : self.nb_nodes_in_riverbed] = np.interp(
                xl_all, np.arange(len(self.section_seq), dtype=float), values_at_sections[:, i_var]
            )
        if self.has_floodplain:
            for (pos_start, pos_end), z in zip(self.casiers_nodes_idx, z_at_casiers):
                values[pos_z, pos_start:pos_end] = z
        return values
