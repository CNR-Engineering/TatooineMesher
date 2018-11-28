"""
Génération d'une surface 2D à partir de la géométrie d'un sous-modèle Crue10
/!\ Pour l'instant seuls les sections profil des branches Saint-Venant sont considérées et sans prise en compte des limites de lits!
"""
from crue10.emh.submodel import SubModel
from crue10.emh.section import SectionProfil
from crue10.utils import CrueError, logger
import numpy as np
import sys

from core.arg_command_line import MyArgParse
from core.base import LigneContrainte, MeshConstructor, SuiteProfilsTravers, ProfilTravers
from core.utils import float_vars


def mesh_crue10_submodel_bathy(args):
    # Read SubModel from xml/shp files
    try:
        model = SubModel(args.infile_etu, args.submodel_name)
        model.read_drso()
        model.read_dptg()
        model.read_shp_traces_sections()
        model.read_shp_branches()
    except FileNotFoundError as e:
        logger.critical(e)
        sys.exit(1)
    except CrueError as e:
        logger.critical(e)
        sys.exit(2)

    #TODO:
    #- merge branches at transitions "SectionIdem/SectionProfil"
    #- generate trace of ProfilIdem

    triangles = {}
    mesh_constr = None

    id_profile = 0
    for i, branche in enumerate(model.iter_on_branches()):
        print("\n\n===== TRAITEMENT DE LA BRANCHE %s =====" % branche.id)
        axe = branche.geom
        profils_travers = SuiteProfilsTravers()
        for section in branche.sections:
            if isinstance(section, SectionProfil):
                coords = list(section.get_coord_3d())
                z_array = np.array([(coord[2],) for coord in coords], dtype=float_vars('Z'))
                profils_travers.add_profile(ProfilTravers(id_profile, [(coord[0], coord[1]) for coord in coords],
                                                          z_array, section.id))
                id_profile += 1

        if len(profils_travers) >= 2:
            profils_travers.compute_dist_proj_axe(axe)
            profils_travers.check_intersections()
            profils_travers.sort_by_dist()
            lignes_contraintes = LigneContrainte.get_lines_from_profils(profils_travers)
            profils_travers.find_and_add_limits(lignes_contraintes, args.dist_max)

            mesh_constr = MeshConstructor(profils_travers, args.pas_trans)
            mesh_constr.build_interp(lignes_contraintes, args.pas_long, args.constant_ech_long)

            if args.outfile_mesh is not None:
                mesh_constr.build_mesh()
                if len(mesh_constr.points) != len(mesh_constr.triangle['vertices']):
                    print("ERROR: Mesh is corrupted... %i vs %i nodes" % (
                          len(mesh_constr.points), len(mesh_constr.triangle['vertices'])))
                else:
                    if not triangles:  # set initial values from first iteration
                        points = mesh_constr.points
                        triangles['triangles'] = mesh_constr.triangle['triangles']
                        triangles['vertices'] = mesh_constr.triangle['vertices']
                    else:  # concatenate with current sub-mesh for next iterations
                        points = np.hstack((points, mesh_constr.points))
                        triangles['triangles'] = np.vstack(
                            (triangles['triangles'], triangles['vertices'].shape[0] + mesh_constr.triangle['triangles']))
                        triangles['vertices'] = np.vstack((triangles['vertices'], mesh_constr.triangle['vertices']))

    mesh_constr.points = points
    mesh_constr.triangle = triangles

    if args.outfile_semis is not None:
        mesh_constr.export_points(args.outfile_semis)
    if args.outfile_mesh is not None:
        mesh_constr.export_mesh(args.outfile_mesh)


parser = MyArgParse(description=__doc__)
# Inputs
parser_infiles = parser.add_argument_group("Données d'entrée du sous-modèle Crue10 à traiter")
parser_infiles.add_argument("infile_etu", help="fichier d'étude Crue10 (etu.xml)")
parser_infiles.add_argument("submodel_name", help="nom du sous-modèle")
# Mesh parameters
parser_mesh = parser.add_argument_group("Paramètres pour la génération du maillage 2D")
parser_mesh.add_argument("--pas_long", type=float, help="pas d'interpolation longitudinal (en m)", default=5)
parser_mesh.add_argument("--pas_trans", type=float, help="pas d'interpolation transversal (en m)", default=3.5)
parser_mesh.add_argument("--dist_max", type=float, help="distance de recherche maxi des 'intersections fictifs' "
                                                        "pour les limites de lits (en m)", default=0.01)
parser_mesh.add_argument("--constant_ech_long",
                         help="méthode de calcul du nombre de profils interpolés entre profils : "
                              "par profil (constant, ie True) ou par lit (variable, ie False)",
                         action='store_true')
# Outputs
parser_outfiles = parser.add_argument_group('Fichiers de sortie')
parser_outfiles.add_argument("--outfile_mesh", help="maillage au format : slf (Telemac), t3s (BlueKenue),"
                                                    " ou xml (LandXML pour ArcGIS)")
parser_outfiles.add_argument("--outfile_semis", help="semis de points au format xyz"
                                                     " (correspond aux noeuds du maillage")


if __name__ == '__main__':
    args = parser.parse_args()
    mesh_crue10_submodel_bathy(args)
