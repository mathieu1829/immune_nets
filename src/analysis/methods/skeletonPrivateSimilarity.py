import numpy as np
import pandas as pd
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.algorithms.simple_beta_distance import simple_beta_distance
from src.creation.distance.levenshtein import levenshteinDistance
import igraph as ig

from src.creation.utils.pathManager import pathManager
from src.creation.io_strategies.test_csv_strategy import *

path = pathManager().testDataPath / "bigTest.csv"


def skeleton_private_similarity(repertoire, minNodeCount=10, minCloneCount=100):
    immuneNet = simple_beta_distance(repertoire = repertoire,distance = levenshteinDistance(group = True),threshold = 2)
    df_net = immuneNet.network
    # print(df_net)
    vertices = np.unique(df_net.to_numpy().flatten())
    net = ig.Graph(df_net.to_numpy())
    net.add_vertices(immuneNet.sampleSize - net.vcount())
    clusters = net.community_fastgreedy()
    clusters = list(clusters.as_clustering(clusters.optimal_count))
    clusters.sort(key = lambda x : len(x))
    private_clusters = clusters[-20:]
    for idx, cluster in enumerate(private_clusters):
        private_clusters[idx] = np.delete(cluster, [ i for i,v in enumerate(cluster) if repertoire.clones.iloc[v]['frequency'] >= minCloneCount] )
    private_clusters = [i for i in private_clusters if len(i) and len(i)>=minNodeCount ]
    cluster_candidates = [e for i in private_clusters for e in i]
    # print(cluster_candidates)

    new_clonotypes = repertoire.clones.iloc[cluster_candidates]
    # print(new_clonotypes)
    new_clonotypes.name = repertoire.clones.name
    new_repertoire = immuneRepertoire(clones=new_clonotypes,sampleIDs={ i:(new_clonotypes["sampleID"].to_numpy() == i).sum() for i in np.unique(new_clonotypes["sampleID"].to_numpy())}) 
    immuneNet = simple_beta_distance(repertoire = new_repertoire,distance = levenshteinDistance(group = True),threshold = 2)
    df_net = immuneNet.network
    # print(df_net)
    vertices = np.unique(df_net.to_numpy().flatten())
    net = ig.Graph(df_net.to_numpy())
    net.add_vertices(immuneNet.sampleSize - net.vcount())
    clusters = net.community_fastgreedy()
    clusters = list(clusters.as_clustering(clusters.optimal_count))
    clusters.sort(key = lambda x : len(x))
    private_clusters = clusters
    # print(private_clusters)
    for idx, cluster in enumerate(private_clusters):
        private_clusters[idx] = np.delete(cluster, [ i for i,v in enumerate(cluster) if repertoire.clones.iloc[v]['frequency'] < minCloneCount] )
    private_clusters = [i for i in private_clusters if len(i) and np.unique(new_clonotypes.iloc[i]["sampleID"]).shape[0] == 1 ]
    return private_clusters

if __name__ == "__main__":
   print(skeleton_private_similarity(test_csv_strategy().input(path)))
