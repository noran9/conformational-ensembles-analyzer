from single_conformation_features import compute_single_conformation_features
from ensemble_features import compute_ensemble_features
from graph import compute_graph
from dendrogram import compute_dendrogram
from pymol_image import generate_image
from mean_rmsd import compute_local_score
from features_images import generate_features_images
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--single", type=str,
                    help="Compute single conformation features for a given .pdb file in the data folder. " +
                         "Add the structure ID as parameter.")
                         
parser.add_argument("-d", "--data", action="store_true",
                    help="Compute single conformation features for all .pdb files in the data folder.")
                    
parser.add_argument("-se", "--single_ensemble", type=str,
                    help="Compute ensemble features for a given .pdb file in the data folder. " +
                         "Add the structure ID as parameter.")
                         
parser.add_argument("-e", "--ensemble", action="store_true",
                    help="Compute ensemble features for all computed single conformation features.")
                    
parser.add_argument("-g", "--graph", type=str,
                    help="Compute similarity graph for a given structure. Add the structure ID as parameter.")
                    
parser.add_argument("-p", type=int, help="Optional argument for the size of the graph. Default value is 0. " +
                    "Use values 1, 2, 3.")

parser.add_argument("-i", "--image", action="store_true",
                    help="Generate PyMOL image from the last computed graph cluster centers.")
                    
parser.add_argument("-dend", "--dendrogram", action="store_true",
                    help="Compute heatmap and dendrogram for similarity among all computed ensembles.")

parser.add_argument("-paired", "--pairwise_dendrogram", nargs='+', default=[],
                    help="Compute heatmap and dendrogram for similarity some ensembles. " +
                         "Add two or more structure IDs as parameters")

parser.add_argument("-l", "--local", action="store_true",
                    help="Compute local score: RMSD variance among all computed ensembles")

parser.add_argument("-fi", "--features_images", action="store_true",
                    help="Compute pairwise comparison plots and summary plots per ensemble")

args = parser.parse_args()

if args.single:
    compute_single_conformation_features(args.single)

if args.data:
    for filename in os.listdir('./data'):
        structure_name = filename.replace('.pdb', '')
        compute_single_conformation_features(structure_name)

if args.single_ensemble:
    compute_ensemble_features(args.single_ensemble)

if args.ensemble:
    for filename in os.listdir('./'):
        if filename.startswith('single-conformations-features-'):
            structure_name = filename.replace('single-conformations-features-', '')
            compute_ensemble_features(structure_name)

if args.graph:
    if args.p:
        compute_graph(args.graph, args.p)
    else:
        compute_graph(args.graph)

if args.dendrogram:
    structure_ids = []
    for filename in os.listdir('./ensemble-features'):
        structure_ids.append(filename)
    compute_dendrogram(structure_ids)

if args.pairwise_dendrogram:
    compute_dendrogram(args.pairwise_dendrogram)

if args.image:
    generate_image()

if args.local:
    compute_local_score()

if args.features_images:
    structure_ids = []
    for filename in os.listdir('./ensemble-features'):
        structure_ids.append(filename)
    generate_features_images(structure_ids)