import argparse


LINE_WIDTH = 80


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


class MyArgParse(argparse.ArgumentParser):
    """
    Custom argparse with a custom formatter_class and a clean help formatter
    Verbose are argument is compulsory

    ### Attributs
    - infile_args
    - mesher_args
    - outfile_args
    """
    def __init__(self, description=None, *args, **kwargs):
        self.infile_args = None
        self.mesher_args = None
        self.outfile_args = None

        kwargs['formatter_class'] = CustomFormatter
        new_description = '_' * LINE_WIDTH + '\n' + description + '_' * LINE_WIDTH + '\n'
        super().__init__(description=new_description, *args, **kwargs)
        self._positionals.title = self._title_group('Positional and compulsory arguments')
        self._optionals.title = self._title_group('Optional arguments')
        self.add_argument("--verbose", "-v", help="increase output verbosity", action="store_true")

    @staticmethod
    def _title_group(label):
        """Decorates group title label"""
        return '~> ' + label

    def add_argument_group(self, name, *args, **kwargs):
        """Add title group decoration"""
        return super().add_argument_group(self._title_group(name), *args, **kwargs)

    def add_common_args(self, project_straight_line=False, constant_long_disc=False):
        self.add_argument("--dist_max", type=float, default=0.01,
                          help="distance de recherche maxi des 'intersections fictifs' "
                               "pour les limites de lits (en m)")

        self.infile_args = self.add_argument_group("Input files")

        self.mesher_args = self.add_argument_group("Mesher and interpolator arguments")
        if project_straight_line:
            self.mesher_args.add_argument("--project_straight_line", action='store_true',
                                          help="project cross-sections along a straight line "
                                               "(linking cross-section bounds)")
        self.mesher_args.add_argument("--long_step", type=float, help="longitudinal space step (in m)", required=True)
        group = self.mesher_args.add_mutually_exclusive_group(required=True)
        group.add_argument("--lat_step", type=float, help="lateral space step (in m)")
        group.add_argument("--nb_pts_lat", type=int, help="number of nodes crosswise")
        if constant_long_disc:
            self.mesher_args.add_argument("--constant_long_disc", action='store_true',
                                          help="method to compute number of intermediate cross-sections "
                                               "per zone (identical between 2 consecutive cross-sections) instead of "
                                               "per bed/submesh (default behaviour)")
        self.mesher_args.add_argument("--interp_constraint_lines", default='LINEAR',
                                      help="interpolation method for X and Y coordinates of constraint lines",
                                      choices=('LINEAR', 'FINITE_DIFF', 'CARDINAL'))
        self.mesher_args.add_argument("--interp_values", default='LINEAR',
                                      help="interpolation method (crosswise for 1D or global 2D) for values",
                                      choices=('LINEAR', 'B-SPLINE', 'AKIMA', 'PCHIP', 'CUBIC_SPLINE',
                                               'BILINEAR', 'BICUBIC', 'BIVARIATE_SPLINE'))

        self.outfile_args = self.add_argument_group("Output files arguments")

    def add_out_mesh_file(self):
        self.outfile_args.add_argument("outfile_mesh", help="Telemac results file")
        self.outfile_args.add_argument("--lang", help="language for standard variables in output file",
                                       default='fr', choices=['en', 'fr'])
