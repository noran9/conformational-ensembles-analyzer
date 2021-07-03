from radius_of_gyration import radius_conformation
from relative_asa import relative_asa
from secondary_structure_project import secondary_structure
from distance_matrix import get_distance_matrix
from Bio.PDB.PDBParser import PDBParser
import os
import numpy as np


def compute_single_conformation_features(structure_name):
    if not os.path.exists('./single-conformations-features-' + structure_name):
        os.makedirs('./single-conformations-features-' + structure_name)

    print('Loading structure...')

    structure = PDBParser(QUIET=True).get_structure(structure_name, "./data/{}.pdb".format(structure_name))

    print('Total number of conformations: ' + str(len(structure)))

    for conformation in range(len(structure)):

        secondary_structure_list = secondary_structure(structure[conformation])

        if len(secondary_structure_list) + 2 != len(structure[conformation]['A']):
          print('Omitting conformation: ' + str(conformation))
          continue

        print('Computing features for conformation ID: ' + str(conformation))

        f = open('./single-conformations-features-' + structure_name + '/conformation-' + str(conformation), 'w')
        f.write('Radius of gyration:\n')
        f.write(str(radius_conformation(structure[conformation]['A'])) + '\n')

        f.write('Relative surface area for each residue: \n')
        asa_list = relative_asa(structure[conformation], structure_name)
        f.write(" ".join(map(str, asa_list)))

        f.write('\nSecondary structure: \n')
        f.write(" ".join(map(str, secondary_structure_list)))

        f.write('\nDistance matrix: \n')
        for line in np.matrix(get_distance_matrix(structure[conformation]['A'])):
            np.savetxt(f, line, fmt='%.2f')

        f.close()