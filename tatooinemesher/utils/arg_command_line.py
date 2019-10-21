import argparse


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


class MyArgParse(argparse.ArgumentParser):
    """
    Custom argparse with a custom formatter_class and a clean help formatter
    Verbose are argument is compulsory
    """
    def __init__(self, *args, **kwargs):
        kwargs['formatter_class'] = CustomFormatter
        super(MyArgParse, self).__init__(*args, **kwargs)
        self.add_argument("--verbose", "-v", help="increase output verbosity", action="store_true")

    def add_common_args(self):
        parser_infiles = self.add_argument_group("Compulsory input files")
        parser_infiles.add_argument("infile_axis", help="fichier d'entrée de l'axe hydraulique (*.shp, *.i2s)")
        parser_infiles.add_argument("infile_cross_sections",
                                    help="fichier d'entrée de profils en travers (*.shp, *.i3s)")
        self.add_argument("--infile_constraint_lines", help="fichier d'entrée de lignes de contrainte (*.shp, *.i2s)")
        self.add_argument("--attr_cross_sections", help="attribut pour identifier les profils en travers")
        self.add_argument("--long_step", type=float, help="longitudinal space step (in m)")
        group = self.add_mutually_exclusive_group(required=True)
        group.add_argument("--lat_step", type=float, help="lateral space step (in m)")
        group.add_argument("--nb_pts_lat", type=int, help="number of nodes crosswise")
        self.add_argument("--dist_max", type=float, help="distance de recherche maxi des 'intersections fictifs' "
                                                         "pour les limites de lits (en m)", default=0.01)
        self.add_argument("--project_straight_line", action='store_true',
                          help="modifie les profils avec une projection plane sur la ligne berge-à-berge")
