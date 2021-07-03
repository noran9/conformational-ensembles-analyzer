import os
import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import ClusterWarning
from warnings import simplefilter


def compute_dendrogram(structure_ids):
    if not os.path.exists('./ensemble-features/'):
        print('Please compute ensemble features for at least two ensembles')
        sys.exit()

    print('Reading ensembles.. ')

    distributions = []
    labels = []

    for filename in os.listdir('ensemble-features/'):
        if filename in structure_ids:
            with open('ensemble-features/' + filename) as f:
                labels.append(filename)

                distance_distribution = []
                next(f)
                next(f)
                next(f)
                next(f)
                next(f)
                next(f)
                next(f)
                next(f)
                next(f)
                while True:
                    line = f.readline().split()
                    if not line:
                        break
                    distance_distribution.append(np.genfromtxt(line))
                distance_distribution = np.asmatrix(distance_distribution)
                f.close()

                new_distance_distribution = []
                for row in distance_distribution:
                    row = np.nan_to_num(row)
                    row = np.asarray(row).reshape(-1)
                    new_distance_distribution.append(row)
                new_distance_distribution = np.asmatrix(new_distance_distribution)
                distributions.append(new_distance_distribution)


    print('Computing global scores..')

    global_score = []
    n = distributions[0][0].shape[1]

    def global_score_calculate(ens1, ens2):
        global_score = np.sqrt(1 / n) * np.linalg.norm(ens1 - ens2)
        return global_score

    for i in range(len(distributions)):
        row = []
        for j in range(len(distributions)):
            if i != j:
                row.append(global_score_calculate(distributions[i], distributions[j]))
            else:
                row.append(0)
        global_score.append(row)
    global_score = np.array(global_score)

    correlations = abs(global_score)

    simplefilter("ignore", ClusterWarning)
    heatmap = sns.clustermap(correlations, annot=True,
                             annot_kws={"size": 7},
                             xticklabels=labels,
                             yticklabels=labels,
                             cmap=sns.cm.rocket_r)

    if not os.path.exists('./images'):
      os.makedirs('./images')

    plt.savefig('./images/clustermap.png')
    print('Heatmap and dendrogram are saved as images/clustermap.png')

