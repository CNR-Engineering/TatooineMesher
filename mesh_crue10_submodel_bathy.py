#!/usr/bin/env python3
"""
Génération d'une surface 2D à partir de la géométrie d'un sous-modèle Crue10
/!\ Pour l'instant seuls les sections profil des branches Saint-Venant sont considérées et sans prise en compte des limites de lits!
"""
from crue10.emh.section import SectionIdem, SectionProfil
from crue10.emh.submodel import SubModel
from crue10.utils import CrueError, logger
import numpy as np
from shapely.affinity import translate
from shapely.geometry import Point
import sys

from core.arg_command_line import MyArgParse
from core.base import LigneContrainte, MeshConstructor, SuiteProfilsTravers, ProfilTravers
from core.utils import float_vars, logger, TatooineException


SHIFT_TRACE_TO_EXTERMITY = True


def mesh_crue10_submodel_bathy(args):
    # Read SubModel from xml/shp files
    try:
        model = SubModel(args.infile_etu, args.submodel_name)
        model.read_drso()
        model.read_dptg()
        model.read_shp_noeuds()
        model.read_shp_traces_sections()
        model.read_shp_branches()
    except FileNotFoundError as e:
        logger.critical(e)
        sys.exit(1)
    except CrueError as e:
        logger.critical(e)
        sys.exit(2)
    logger.info(model)

    #TODO:
    # - merge branches at transitions "SectionIdem/SectionProfil"
    # - generate trace of ProfilIdem

    if SHIFT_TRACE_TO_EXTERMITY:
        # Shift first trace at branch origin (usual biais is about 15.039 m)
        for i, branche in enumerate(model.iter_on_branches()):
            for pos in (0, -1):
                section = branche.sections[pos]
                if isinstance(section, SectionProfil):
                    if pos == 0:
                        node = branche.noeud_amont.geom
                    else:
                        node = branche.noeud_aval.geom
                    section_point = section.geom_trace.intersection(branche.geom)
                    if isinstance(section_point, Point):
                        dx = node.x - section_point.x
                        dy = node.y - section_point.y
                        model.branches[branche.id].sections[pos].set_trace(
                            translate(section.geom_trace, xoff=dx, yoff=dy))

    for branche in model.iter_on_branches():
        for j, section in enumerate(branche.sections):
            if isinstance(section, SectionIdem):
                if j == 0 or j == len(branche.sections) - 1:
                    # Determine node name at junction
                    nom_noeud = branche.noeud_amont.id if j == 0 else branche.noeud_aval.id

                    # Check if any adjacent branchs has this section
                    branches = model.connected_branches(nom_noeud)
                    branches.remove(branche)
                    if branches:
                        for br in branches:
                            section_pos = 0 if br.noeud_amont.id == nom_noeud else -1
                            section_at_node = br.sections[section_pos]
                            if section_at_node is section.section_ori:
                                branche.sections[j] = section.set_as_profil(section_at_node, section.dz)
                            else:
                                pass  # TODO: build trace
                else:
                    pass  #TODO: build trace

    triangles = {}
    mesh_constr = None

    id_profile = 0
    for i, branche in enumerate(model.iter_on_branches()):
        logger.info("===== TRAITEMENT DE LA BRANCHE %s =====" % branche.id)
        axe = branche.geom
        try:
            profils_travers = SuiteProfilsTravers()
            for section in branche.sections:
                if isinstance(section, SectionProfil):
                    coords = list(section.get_coord_3d())
                    z_array = np.array([(coord[2],) for coord in coords], dtype=float_vars('Z'))
                    profils_travers.add_profile(ProfilTravers(id_profile, [(coord[0], coord[1]) for coord in coords],
                                                              z_array, section.id))
                    id_profile += 1

            if len(profils_travers) >= 2:
                profils_travers.compute_dist_proj_axe(axe, args.dist_max)
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
            else:
                logger.info("Branche ignorée")
            logger.info("\n\n")
        except TatooineException as e:
            logger.critical(e.message)

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
