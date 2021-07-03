import matplotlib.pyplot as plt
import os
import numpy as np
import seaborn as sns


def compute_local_score():
    rmsd = []

    for filename in os.listdir('./ensemble-features/'):
        with open('./ensemble-features/' + filename, 'r') as f:

            print('Reading conformations {}... '.format(filename))
            next(f)
            next(f)
            next(f)
            next(f)
            next(f)
            next(f)
            next(f)

            mrmsd = f.readline().strip().split(' ')
            rmsd.append([float(x) for x in mrmsd])

            f.close()

    # Compute variance per residue
    rmsd_line = np.var(rmsd, axis=0)

    # Normalize array
    rmsd_line = (rmsd_line - np.min(rmsd_line)) / (np.max(rmsd_line) - np.min(rmsd_line))

    # Plot the variance of RMSD among ensembles
    sns.set_style("darkgrid")
    plt.figure(figsize=(50, 15))
    plt.plot(rmsd_line, c='orange', ls='--', lw=7)
    plt.xlabel("Residue ID", size=30)
    plt.ylabel("Normalized RMSD Variance", size=30)
    plt.xticks(range(len(rmsd_line)), fontsize=20)
    plt.yticks(fontsize=20)
    plt.title("Variance of local RMSD metric per residue", size=40)

    if not os.path.exists('./images'):
      os.makedirs('./images')

    plt.savefig('images/local-rmsd-score.png', bbox_inches='tight')
    plt.close()

    print("Local score plot generated as /images/local-rmsd-score.png.")

