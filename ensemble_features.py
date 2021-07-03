from Bio.PDB.PDBParser import PDBParser
from RMSD import get_average_rmsd
import os
import sys
import numpy as np


def compute_ensemble_features(structure_name):
    secondary_structure_class = ['E', 'P', 'H', 'L']

    if not os.path.exists('./ensemble-features'):
        os.makedirs('./ensemble-features')

    if not os.path.exists('./single-conformations-features-' + structure_name):
        print('Please compute single conformation features.')
        sys.exit()

    ensemble_radius, ensemble_ss, ensemble_msa, ensemble_mat = [], [], [], []

    print('Reading conformations.. ')

    for filename in os.listdir('single-conformations-features-' + structure_name + '/'):
        with open('./single-conformations-features-' + structure_name + '/' + filename, 'r') as f:
            next(f)
            ensemble_radius.append(float(f.readline().strip()))

            next(f)
            asa = f.readline().split()
            ensemble_msa.append([float(x) for x in asa])

            next(f)
            ss = f.readline().split()
            ensemble_ss.append([secondary_structure_class.index(x) for x in ss])

            next(f)
            model_matrix = []
            while True:
                line = f.readline().split()
                if not line:
                    break
                model_matrix.append(np.genfromtxt(line))
            model_matrix = np.asmatrix(model_matrix)
            ensemble_mat.append(model_matrix)

            f.close()

    print('Writing to file.. ')

    f = open('./ensemble-features/' + structure_name, 'w')

    f.write('Radius of gyration:\n')

    for rad in ensemble_radius:
        f.write(str(round(rad, 2)) + " ")

    f.write('\nMedian solvent accessibility: \n')

    for s in np.average(ensemble_msa, axis=0):
        f.write(str(round(s, 2)) + " ")

    f.write('\nSecondary structure entropy: \n')
    ensemble_ss = np.average(ensemble_ss, axis=0)
    ensemble_ss = [secondary_structure_class[round(x)] for x in ensemble_ss]
    f.write(" ".join(map(str, ensemble_ss)))

    f.write('\nMedian RMSD: \n')

    structure = PDBParser(QUIET=True).get_structure(structure_name, "./data/{}.pdb".format(structure_name))
    for rmsd in get_average_rmsd(structure):
        f.write(str(round(rmsd, 2)) + ' ')

    f.write('\nMean distance matrix: \n')
    np.savetxt(f, np.average(ensemble_mat, axis=0), fmt='%.2f')

    f.write('\nStandard deviation of distance matrix: \n')
    np.savetxt(f, np.std(ensemble_mat, axis=0), fmt='%.2f')

    f.close()
    print('Ensemble features are generated for ID: ' + structure_name)
