import numpy as np
import pandas as pd
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.algorithms.simple_beta_distance import simple_beta_distance
from src.creation.distance.levenshtein import levenshteinDistance
from src.analysis.methods.publicSimilarity import publicSimilarity
import igraph as ig
import itertools

from src.creation.utils.pathManager import pathManager
from src.creation.io_strategies.test_csv_strategy import *

path = pathManager().testDataPath / "bigTest.csv"


def privateSimilarity(repertoire, frequencyCutoff = None, top_k = 20, absoulutePublic=False, minCoverage = 2, minClusterSize=2, minPrivateClusterSize=1):
    public_clusters = publicSimilarity(repertoire=repertoire, 
                                       frequencyCutoff=frequencyCutoff,
                                       top_k=top_k,absoulutePublic=absoulutePublic,
                                       minCoverage=minCoverage,
                                       minClusterSize=minClusterSize)
    
    public_clusters = set(list(itertools.chain(*public_clusters)))

    prepared_clones = repertoire.clones.dropna(subset = ["tcra_aa", "tcrb_aa"]) 
    prepared_clones.name = repertoire.clones.name
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
    private_clusters = clusters[-top_k:]
    # print(public_clusters)


    # assigning indices from original df
    private_clusters = [ prepared_clones.iloc[cluster].index.tolist() for cluster in private_clusters]


    # collecting clusters that have no overlap with public and have big enough size
    private_clusters = [ cluster for cluster in private_clusters if not set(cluster) & public_clusters and len(cluster) >= minPrivateClusterSize and len(np.unique(prepared_clones.loc[cluster]["sampleID"])) == 1]

    clusterSampleIDs = np.array([ prepared_clones.loc[cluster[0]]["sampleID"] for cluster in private_clusters])
    private_clusters = np.array(private_clusters)

    #separating clusters by sampleIDs
    private_clusters = [ private_clusters[clusterSampleIDs == sampleID].tolist() for sampleID in np.unique(prepared_clones["sampleID"])]
    return private_clusters


if __name__ == "__main__":
   print(privateSimilarity(test_csv_strategy().input(path)))
