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
        parser_infiles = self.add_argument_group("Fichiers d'entrée obligatoires")
        parser_infiles.add_argument("infile_axe", help="fichier d'entrée de l'axe hydraulique (i2s ou shp)")
        parser_infiles.add_argument("infile_profils_travers", help="fichier d'entrée de profils en travers (i3s ou shp)")
        self.add_argument("--infile_lignes_contraintes", help="fichier d'entrée de lignes de contrainte (i2s ou shp)")
        self.add_argument("--attr_profils_travers", help="attribut pour identifier les profils en travers")
        self.add_argument("--pas_long", type=float, help="pas d'interpolation longitudinal (en m)", default=5)
        self.add_argument("--pas_trans", type=float, help="pas d'interpolation transversal (en m)", default=3.5)
        self.add_argument("--dist_max", type=float, help="distance de recherche maxi des 'intersections fictifs' "
                                                         "pour les limites de lits (en m)", default=0.01)
