import numpy as np
import pandas as pd
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.algorithms.simple_beta_distance import simple_beta_distance
from src.creation.distance.levenshtein import levenshteinDistance
import igraph as ig

from src.creation.utils.pathManager import pathManager
from src.creation.io_strategies.test_csv_strategy import *

path = pathManager().testDataPath / "bigTest.csv"


def publicSimilarity(repertoire, frequencyCutoff = None, top_k = 20, absoulutePublic=False, minCoverage = 2, minClusterSize=2):
    prepared_clones = repertoire.clones.dropna(subset = ["tcra_aa", "tcrb_aa"]) 
    prepared_clones.name = repertoire.clones.name

    # If specified take into account only clones above frequency cutoff
    if frequencyCutoff is not None:
        prepared_clones[prepared_clones["frequency"] > frequencyCutoff]

    preparedRepertoire = immuneRepertoire(clones=prepared_clones)
    # print(prepared_clones)

    
    immuneNet = simple_beta_distance(repertoire = preparedRepertoire,distance = levenshteinDistance(group = True),threshold = 2)
    df_net = immuneNet.network

    #performing clustering
    vertices = np.unique(df_net.to_numpy().flatten())
    net = ig.Graph(df_net.to_numpy())
    net.add_vertices(immuneNet.sampleSize - net.vcount())
    clusters = net.community_fastgreedy()
    clusters = list(clusters.as_clustering(clusters.optimal_count))

    # print(prepared_clones.iloc[clusters[0][0]]["frequency"])
    #filtering out k=20 largerst clusters
    clusters.sort(key = lambda x : sum([ prepared_clones.iloc[i]["frequency"] for i in x]))
    # clusters.sort(key = lambda x : len(x))
    # for i in clusters:
    #     print(i)
    public_clusters = clusters[-top_k:]
    # print(public_clusters)

    # assigning indices from original df
    public_clusters = [ prepared_clones.iloc[cluster].index.tolist() for cluster in public_clusters]

    #fitering for clusters with size bigger than min size and by how many record samples are there 
    public_clusters = [ cluster for cluster in public_clusters if (len(cluster) >= minClusterSize and len(np.unique(prepared_clones.loc[cluster]["sampleID"])) >= (minCoverage if not absoulutePublic else len(np.unique(prepared_clones["sampleID"])))) ]

    return public_clusters
