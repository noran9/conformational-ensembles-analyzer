import math
import numpy as np
from Bio.PDB.PDBParser import PDBParser


# Radius of gyration
def radius_conformation(chain):

    # Heavy atoms coordinates
    coord = list()
    for atom in chain.get_atoms():
        if atom.get_name()[0] in ['C', 'O', 'N', 'S']:
            coord.append(atom.get_coord())
    coord = np.array(coord)  # N X 3

    barycenter = np.sum(coord, axis=0) / coord.shape[0]  # center of mass is more correct

    # Calculate distance of each atom from the barycenter
    dist = coord - barycenter
    dist = dist * dist
    dist = np.sqrt(np.sum(dist, axis=1))

    return round(math.sqrt(np.sum(dist * dist) / len(coord)), 3)


def radius_ensemble(structure_id):
    # Load the structure
    structure = PDBParser(QUIET=True).get_structure(structure_id, "./data/{}.pdb".format(structure_id))

    radius = []
    for i in range(len(structure)):
        radius.append(radius_conformation(structure[i]['A']))

    return radius
