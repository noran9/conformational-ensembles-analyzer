import matplotlib.pyplot as plt
import os
import sys
import numpy as np
import seaborn as sns


def generate_features_images(structure_names):
    secondary_structure_class = ['E', 'P', 'H', 'L']
    first_residue = 11

    if not os.path.exists('./local-score-images'):
        os.makedirs('./local-score-images')

    median_rmsd, ensemble_radius, ensemble_ss_entropy, ensemble_msa, ensemble_mat, ensemble_mat_ds = [], [], [], [], [], []

    for structure in structure_names:
        if not os.path.exists('./ensemble-features/' + structure):
            print('Please compute ensemble features.')
            sys.exit()

        print('Reading conformations {}... '.format(structure))

        with open('./ensemble-features/' + structure, 'r') as f:
            next(f)
            er = f.readline().strip().split(' ')
            ensemble_radius.append(np.genfromtxt(er))

            next(f)
            asa = f.readline().split()
            ensemble_msa.append([float(x) for x in asa])

            next(f)
            ss = f.readline().split()
            ensemble_ss_entropy.append([secondary_structure_class.index(x) for x in ss])

            next(f)
            mrmsd = f.readline().strip().split(' ')
            median_rmsd.append([float(x) for x in mrmsd])


            next(f)
            model_matrix = []
            while True:
                line = f.readline().split()
                if not line:
                    break
                model_matrix.append(np.genfromtxt(line))
            model_matrix = np.asmatrix(model_matrix)
            ensemble_mat.append(model_matrix)

            next(f)
            matrix = []
            while True:
                line = f.readline().split()
                if not line:
                    break
                matrix.append(np.genfromtxt(line))
            matrix = np.asmatrix(matrix)
            ensemble_mat_ds.append(matrix)
            f.close()

    new_matrixfin = []
    for i in range(len(structure_names)):
        matrixfin = np.triu(ensemble_mat[i]) + np.tril(ensemble_mat_ds[i])
        new_matrix = []
        for row in matrixfin:
            row = np.nan_to_num(row)
            row = np.asarray(row).reshape(-1)
            new_matrix.append(row)
        new_matrix = np.asmatrix(new_matrix)
        new_matrixfin.append(new_matrix)


    ############################ PLOT ALL ENSEMBLE FEATURES IN A IMAGE ############################
    for i in range(len(structure_names)):
        x_pos = np.arange(first_residue, len(median_rmsd[i]) + first_residue)

        plt.figure(figsize=(25, 20))
        plt.plot()

        rmsd = plt.subplot(2, 2, 1)
        msa = plt.subplot(2, 2, 2)
        rg = plt.subplot(2, 2, 3)
        dm = plt.subplot(2, 2, 4)

        rmsd.plot(x_pos, median_rmsd[i])
        rmsd.set_xlabel('Residues', fontsize=25)
        rmsd.set_ylabel('RMSD', fontsize=25)
        rmsd.tick_params(axis='both', labelsize=15)

        msa.plot(x_pos, ensemble_msa[i])
        msa.set_xlabel('Residues', fontsize=25)
        msa.set_ylabel('MSA', fontsize=25)
        msa.tick_params(axis='both', labelsize=15)

        sns.histplot(ensemble_radius[i], ax=rg)
        rg.set_xlabel('Radius of gyration', fontsize=25)
        rg.set_ylabel('Count', fontsize=25)
        rg.tick_params(axis='both', labelsize=15)

        dm.imshow(new_matrixfin[i])
        dm.tick_params(axis='both', labelsize=15)

        plt.suptitle(t='Features of ensemble {}'.format(structure_names[i]), fontsize=35)

        plt.savefig('./local-score-images/Features_{}.png'.format(structure_names[i]), bbox_inches='tight')

        plt.close()

        print("Files generated in /local-score-images.")
