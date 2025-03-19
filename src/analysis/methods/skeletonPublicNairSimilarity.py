import numpy as np
import pandas as pd
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.algorithms.simple_beta_distance import simple_beta_distance
from src.creation.distance.levenshtein import levenshteinDistance
import igraph as ig

from src.creation.utils.pathManager import pathManager
from src.creation.io_strategies.test_csv_strategy import *

path = pathManager().testDataPath / "bigTest.csv"


def skeletonPublicNairSimilarity(repertoire, top_k = 20, absoulutePublic=False, minCoverage = 2, minClusterSize=2):
    prepared_clones = repertoire.clones.dropna(subset = ["tcra_aa", "tcrb_aa"]) 
    # print(prepared_clones)

    dictIdx = {}
    skeleton_clones = pd.DataFrame()
    # Creating network and cluster analysis for each sample
    for sampleID in np.unique(prepared_clones["sampleID"]):
        # print(f"sampleID: {sampleID}")
        #creating network for selected sample
        sampleClones = prepared_clones.iloc[ prepared_clones["sampleID"].to_numpy() == sampleID]
        sampleClones.name = repertoire.clones.name
        sampleRepertoire = immuneRepertoire(clones=sampleClones)
        immuneNet = simple_beta_distance(repertoire = sampleRepertoire,distance = levenshteinDistance(group = True),threshold = 2)
        df_net = immuneNet.network

        #performing clustering
        vertices = np.unique(df_net.to_numpy().flatten())
        net = ig.Graph(df_net.to_numpy())
        net.add_vertices(immuneNet.sampleSize - net.vcount())
        clusters = net.community_fastgreedy()
        clusters = list(clusters.as_clustering(clusters.optimal_count))

        # print(prepared_clones.iloc[clusters[0][0]]["frequency"])
        #filtering out k=20 largerst clusters
        clusters.sort(key = lambda x : sum([ sampleClones.iloc[i]["frequency"] for i in x]))
        # clusters.sort(key = lambda x : len(x))
        # for i in clusters:
        #     print(i)
        public_clusters = clusters[-top_k:]
        # print(public_clusters)

        #picking representatives aka skeleton clones
        clusterRepresentativeClones = [ cluster[sampleClones.iloc[cluster]['frequency'].to_numpy().argmax()] for cluster in public_clusters]
        sample_skeleton_clones = sampleClones.iloc[clusterRepresentativeClones]

        #assigning indexes for cluster a representative was derived from
        sample_skeleton_clones["clusterIdx"] = np.arange(20)

        public_clusters = [ sampleClones.iloc[cluster].index.tolist() for cluster in public_clusters]
        dictIdx[sampleID] = public_clusters

        # print(sample_skeleton_clones.iloc[:,1:])
        skeleton_clones = pd.concat([skeleton_clones,sample_skeleton_clones.iloc[:,1:]])

    # Creating a network from skeleton clones
    # print(skeleton_clones)
    skeleton_clones.name = repertoire.clones.name
    skeleton_repertoire = immuneRepertoire(clones=skeleton_clones)
    immuneNet = simple_beta_distance(repertoire = skeleton_repertoire,distance = levenshteinDistance(group = True),threshold = 2)
    df_net = immuneNet.network

    #performing clustering on skeleton clones
    vertices = np.unique(df_net.to_numpy().flatten())
    net = ig.Graph(df_net.to_numpy())
    net.add_vertices(immuneNet.sampleSize - net.vcount())
    clusters = net.community_fastgreedy()
    clusters = list(clusters.as_clustering(clusters.optimal_count))

    
    clusters.sort(key = lambda x : sum([ skeleton_clones.iloc[i]["frequency"] for i in x]))

    # assigning the indexes from original dfs to clusters
    clusters = [ skeleton_clones.iloc[cluster].index.tolist() for cluster in clusters]

    # discarding clusters with records from only one sample
    clusters = [ cluster for cluster in clusters if (len(cluster) >= minClusterSize and len(np.unique(skeleton_clones.loc[cluster]["sampleID"])) >= (minCoverage if not absoulutePublic else len(np.unique(prepared_clones["sampleID"]))) ) ]
    
    # for i,cluster in enumerate(clusters):
    #     print(f"cluster no. {i}:") 
    #     for dfIdx in cluster:
    #         print(f"{dfIdx}: {skeleton_clones.loc[dfIdx]["sampleID"]}")

    # expanding the skeleton
    for cluster in clusters:
        for idx in range(len(cluster)):
            dfIdx = cluster[idx]
            # print(f"dfIdx: {dfIdx}")
            initialCluster = dictIdx[prepared_clones.loc[dfIdx]["sampleID"]][skeleton_clones.loc[dfIdx]["clusterIdx"]]
            # print(f"initialCluster: {initialCluster}")
            initialCluster.remove(dfIdx)
            # print(f"initialCluster post: {initialCluster}")
            cluster.extend(initialCluster)

    #public clusters ready
    return clusters
    # for i,cluster in enumerate(clusters):
    #     print(f"cluster no. {i}:") 
    #     for dfIdx2 in cluster:
    #         print(f"{dfIdx2}: {prepared_clones.loc[dfIdx2]["sampleID"]}")

if __name__ == "__main__":
   print(skeletonPublicNairSimilarity(test_csv_strategy().input(path),absoulutePublic=True))
