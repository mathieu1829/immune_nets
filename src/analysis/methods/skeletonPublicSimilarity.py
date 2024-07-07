import numpy as np
import pandas as pd
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.algorithms.simple_beta_distance import simple_beta_distance
from src.creation.distance.levenshtein import levenshteinDistance
import igraph as ig

from src.creation.utils.pathManager import pathManager
from src.creation.io_strategies.test_csv_strategy import *

path = pathManager().testDataPath / "bigTest.csv"


def skeleton_similarity(repertoire, minNodeCount=1, minCloneCount=1, absoulute_public=False, minimum_coverage = 2):
    immuneNet = simple_beta_distance(repertoire = repertoire,distance = levenshteinDistance(group = True),threshold = 2)
    df_net = immuneNet.network
    # print(df_net)
    vertices = np.unique(df_net.to_numpy().flatten())
    net = ig.Graph(df_net.to_numpy())
    net.add_vertices(immuneNet.sampleSize - net.vcount())
    clusters = net.community_fastgreedy()
    clusters = list(clusters.as_clustering(clusters.optimal_count))
    clusters.sort(key = lambda x : len(x))
    public_clusters = clusters[-20:]
    # print(public_clusters)
    prepared_clones = repertoire.clones.dropna(subset = ["tcra_aa", "tcrb_aa"])
    for idx, cluster in enumerate(public_clusters):
        public_clusters[idx] = np.delete(cluster, [ i for i,v in enumerate(cluster) if prepared_clones.iloc[v]['frequency'] < minCloneCount] )
    public_clusters = [i for i in public_clusters if len(i) and len(i)>=minNodeCount ]
    # print(public_clusters)
    if absoulute_public == True :
        num_of_samples = np.unique(prepared_clones["sampleID"].to_numpy()).shape[0]
        public_clusters = [i for i in public_clusters if np.unique(prepared_clones.iloc[i]["sampleID"].to_numpy()).shape[0] == num_of_samples]
    else:
        public_clusters = [i for i in public_clusters if np.unique(prepared_clones.iloc[i]["sampleID"].to_numpy()).shape[0] >= minimum_coverage]
    
    return public_clusters

if __name__ == "__main__":
   print(skeleton_similarity(test_csv_strategy().input(path), absoulute_public = True))
