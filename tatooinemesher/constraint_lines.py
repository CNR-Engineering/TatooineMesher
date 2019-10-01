import numpy as np
from pyteltools.geom import BlueKenue as bk, Shapefile as shp
import shapefile
from shapely.geometry import LineString

from tatooinemesher.interp.cubic_hermite_spline import CubicHermiteSpline
from tatooinemesher.utils import float_vars


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
