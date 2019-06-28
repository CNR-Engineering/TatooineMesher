#!/usr/bin/env python3
"""
Génération d'une surface 2D à partir de la géométrie d'un sous-modèle Crue10
Variables :
* FOND
* FROTTEMENT
* IS_LIT_ACTIF

Branches: seules les branches avec des SectionProfil et/ou SectionIdem sont traitées :
- EMHBrancheSaintVenant (type 20)
- EMHBrancheSeuilTransversal (type 4)
- EMHBrancheStrickler (type 6)
- EMHBarrageFilEau (type 15)
L'option `--types_branches` permet de sélectionner les types de branches à traiter.

Particularités :
- SectionProfil: les sections situées aux extrémités (amont/aval) des branches sont replacés pour correspondre
    géographiquement au début/fin de la branche
- SectionIdem:
    - l'abscisse curviligne xp est normalisée (par la géométrie)
    - la traceProfil est copiée (si la SectionProfil originale est au même endroit) ou regénérée (perpendiculaire à la branche)
- Les branches non actives sont également traitées
"""
import numpy as np
import sys
from time import perf_counter

from crue10.emh.branche import Branche
from crue10.emh.section import SectionProfil
from crue10.emh.submodel import SubModel
from crue10.utils import CrueError, logger

from core.arg_command_line import MyArgParse
from core.base import LigneContrainte, MeshConstructor, SuiteProfilsTravers, ProfilTravers
from core.utils import logger, set_logger_level, TatooineException


def mesh_crue10_submodel(args):
    set_logger_level(args.verbose)
    t1 = perf_counter()

    # Read SubModel from xml/shp files
    try:
        submodel = SubModel(args.infile_etu, args.submodel_name)
        submodel.read_all()
        submodel.remove_sectioninterpolee()
        submodel.normalize_geometry()
        # submodel.write_shp_active_trace('export_tracesSections.shp')  # DEBUG
    except FileNotFoundError as e:
        logger.critical(e)
        sys.exit(1)
    except CrueError as e:
        logger.critical(e)
        sys.exit(2)
    logger.info(submodel)

    triangles = {}
    mesh_constr = None
    points = None
    suite = []

    id_profile = 0
    for i, branche in enumerate(submodel.iter_on_branches()):
        logger.info("===== TRAITEMENT DE LA BRANCHE %s =====" % branche.id)
        axe = branche.geom
        try:
            profils_travers = SuiteProfilsTravers()
            for section in branche.sections:
                if isinstance(section, SectionProfil):
                    coords = list(section.get_coord(add_z=True))
                    profile = ProfilTravers(section.id, [(coord[0], coord[1]) for coord in coords], 'Section')
                    z = np.array([coord[2] for coord in coords])
                    is_bed_active = section.get_is_bed_active_array()
                    mean_strickler = section.get_friction_coeff_array()
                    profile.coord.values = np.core.records.fromarrays(
                        np.column_stack((z, is_bed_active, mean_strickler)).T,
                        names=['B', 'IS BED ACTIVE', 'W']
                    )
                    profils_travers.add_profile(profile)
                    id_profile += 1

            if len(profils_travers) >= 2:
                profils_travers.compute_dist_proj_axe(axe, args.dist_max)
                profils_travers.check_intersections()
                # profils_travers.sort_by_dist() is useless because profiles are already sorted
                lignes_contraintes = LigneContrainte.get_lines_from_profils(profils_travers)
                profils_travers.find_and_add_limits(lignes_contraintes, args.dist_max)

                mesh_constr = MeshConstructor(profils_travers, args.pas_trans)
                mesh_constr.build_interp(lignes_contraintes, args.pas_long, args.constant_ech_long)

                if args.outfile_mesh is not None:
                    mesh_constr.build_mesh()
                    suite += mesh_constr.profils_travers
                    if not triangles:  # set initial values from first iteration
                        points = mesh_constr.points
                        triangles['triangles'] = mesh_constr.triangle['triangles']
                        triangles['vertices'] = mesh_constr.triangle['vertices']
                    else:  # concatenate with current sub-mesh for next iterations
                        last_zone = points['zone'][-1]
                        mesh_constr.points['zone'] += last_zone + 2
                        points = np.hstack((points, mesh_constr.points))
                        triangles['triangles'] = np.vstack(
                            (triangles['triangles'],
                             triangles['vertices'].shape[0] + mesh_constr.triangle['triangles']))
                        triangles['vertices'] = np.vstack((triangles['vertices'], mesh_constr.triangle['vertices']))
            else:
                logger.info("Branche ignorée par manque de sections")
        except TatooineException as e:
            logger.critical("/!\ Branche ignorée à cause d'une erreur bloquante :")
            logger.critical(e.message)
        logger.info("\n")

    mesh_constr.points = points
    mesh_constr.profils_travers = suite
    mesh_constr.triangle = triangles

    logger.info(mesh_constr.summary())  # General information about the merged mesh
    if args.outfile_semis is not None:
        mesh_constr.export_points(args.outfile_semis)
    if args.outfile_mesh is not None:
        mesh_constr.export_mesh(args.outfile_mesh)

    t2 = perf_counter()
    logger.info("=> le temps d'execution est de : {}s".format(t2-t1))


parser = MyArgParse(description=__doc__)
# Inputs
parser_infiles = parser.add_argument_group("Données d'entrée du sous-modèle Crue10 à traiter")
parser_infiles.add_argument("infile_etu", help="fichier d'étude Crue10 (etu.xml)")
parser_infiles.add_argument("submodel_name", help="nom du sous-modèle")
help_branche_types = "Type des branches à considérer (liste d'entiers) : "
for i, id_type in enumerate(Branche.TYPES_WITH_GEOM):
    if i != 0: help_branche_types += ', '
    help_branche_types += '%i = %s' % (id_type, Branche.TYPES[id_type])
parser_infiles.add_argument("--types_branches", help=help_branche_types, type=int, nargs='+',
                            choices=Branche.TYPES_WITH_GEOM, default=[20])
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
    mesh_crue10_submodel(args)
