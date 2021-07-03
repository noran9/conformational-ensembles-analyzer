import networkx as nx
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import os
import sys


def compute_graph(structure_name, p=0):
    if not os.path.exists('./single-conformations-features-' + structure_name):
        print('Please compute single conformation features.')
        sys.exit()

    secondary_structure_class = ['E', 'P', 'H', 'L']

    # Calculate feature lists
    print('Loading features.. ')

    # Load each feature file as row in the features matrix
    features = []
    for filename in os.listdir('single-conformations-features-' + structure_name + '/'):
        with open('./single-conformations-features-' + structure_name + '/' + filename, 'r') as f:

            # Initialize conformation list
            conformation_features = []

            # # Read radius of gyration
            next(f)
            conformation_features.append(float(f.readline().strip()))

            next(f)
            next(f)

            # Read secondary structure
            next(f)
            ss = f.readline().split()
            for x in ss:
                # Map secondary structure characters to numbers
                # Keeping secondary structure values as integer influences the graph partitioning too much
                # The values are normalized in the range of values for secondary structure
                conformation_features.append(secondary_structure_class.index(x) / 3)

            # Flatten the array for this file and add to features matrix
            conformation_features = np.asarray(conformation_features)
            conformation_features.flatten()
            features.append(conformation_features)

    print('Calculate optimal K.. ')

    inertia = []
    centroids = []

    # Number of K values to be checked
    maximum_clusters = 7
    #
    # Run K-Means algorithm for K from 1 to 7
    for k in range(1, maximum_clusters + 1):
        print('K = ' + str(k))
        kmeans = KMeans(n_clusters=k).fit(features)
        inertia.append(kmeans.inertia_)
        centroids.append(kmeans.cluster_centers_)

    if not os.path.exists('./images'):
      os.makedirs('./images')

    # Plot inertia
    plt.plot(inertia)
    plt.xticks(range(maximum_clusters), labels=range(1, maximum_clusters + 1))
    plt.savefig('./images/k-means-inertia.png')
    plt.close()
    print('Inertia plot for different K values is saved as images/k-means-inertia.png')

    # Locate function knee for optimal K
    point_a = np.array([0, inertia[0]])
    point_b = np.array([maximum_clusters - 1, inertia[maximum_clusters - 1]])

    max_distance, optimal_k = 0, 0
    for k in range(1, maximum_clusters + 1):
        point_c = np.array([k, inertia[k-1]])

        d = abs(np.cross(point_b - point_a, point_c - point_a) / np.linalg.norm(point_b - point_a))
        if d > max_distance:
            max_distance = d
            optimal_k = k

    # Return optimal clusters
    print('Optimal number of clusters: ' + str(optimal_k))

    # Conformation closest to the cluster centers
    closest, _ = pairwise_distances_argmin_min(centroids[optimal_k - 1], features)
    print('Representative conformations closest to cluster centers: ')
    print(closest)

    # Generate graph nodes starting with the conformations closest to the centroids
    nodes = list(closest)
    vectors, removed = [], []

    remaining_features = list(features)
    # Iteratively add similar nodes to the graph
    for _ in range(p):
        for c in nodes:
            if c not in removed and c < len(remaining_features):
                vectors.append(remaining_features[c])
                remaining_features.pop(c)
                removed.append(c)
        neighbors, _ = pairwise_distances_argmin_min(vectors, remaining_features)
        nodes = np.append(nodes, neighbors)

    # Weights calculation:
    G = nx.Graph()
    for node_i in range(len(nodes)):
        for node_j in range(len(nodes)):
            # Weights are computed as difference between feature lists between two nodes
            dist = sum(abs(np.array(features[node_i]) - np.array(features[node_j])))

            # Avoid division with 0 and use the reciprocal distance to keep similar nodes closer
            weight = dist
            # Add edge to graph and avoid self loops
            if node_i != node_j:
                G.add_edge(node_i, node_j, weight=weight)

    # Color the nodes according to cluster
    kmeans = KMeans(n_clusters=optimal_k).fit(features)
    colors = []
    for node in nodes:
        colors.append(kmeans.predict([features[node]]))

    # Draw the graph
    # Possible layout functions : nx.spectral_layout; nx.spring_layout (unweighted);
    # nx.shell_layout (unweighted); nx.kamada_kawai_layout
    pos = nx.kamada_kawai_layout(G)

    # Color the nodes according to their partition
    cmap = cm.get_cmap('viridis', len(colors))

    nx.draw(G, pos, node_size=60, node_color=colors, cmap=cmap)

    nx.draw_networkx_edges(G, pos, alpha=0.2)

    plt.savefig('./images/graph-' + structure_name + '.png')

    print('Generated graph image: graph-' + structure_name + '.png in /images.')

    # Save cluster centers used to generate a PyMOL image in a file
    f = open('./cluster-centers', 'w')
    f.write(structure_name + '\n')
    f.write(str(optimal_k) + '\n')
    for c in closest:
        f.write(str(c) + " ")
    f.close()
