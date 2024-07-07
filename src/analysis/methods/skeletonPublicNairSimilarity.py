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
    prepared_clones = repertoire.clones.dropna(subset = ["tcra_aa", "tcrb_aa"]) 
    df_list = []
    dictIdx = {}
    for sampleID in np.unique(prepared_clones["sampleID"]):
        sampleClones = prepared_clones.iloc[ prepared_clones["sampleID"].to_numpy() == sampleID]
        sampleClones.name = repertoire.clones.name
        sampleRepertoire = immuneRepertoire(clones=sampleClones,sampleIDs={ i:(sampleClones["sampleID"].to_numpy() == i).sum() for i in np.unique(sampleClones["sampleID"].to_numpy())})
        immuneNet = simple_beta_distance(repertoire = sampleRepertoire,distance = levenshteinDistance(group = True),threshold = 2)
        df_net = immuneNet.network
        vertices = np.unique(df_net.to_numpy().flatten())
        net = ig.Graph(df_net.to_numpy())
        net.add_vertices(immuneNet.sampleSize - net.vcount())
        clusters = net.community_fastgreedy()
        clusters = list(clusters.as_clustering(clusters.optimal_count))
        clusters.sort(key = lambda x : len(x))
        public_clusters = clusters[-20:]
        # print(public_clusters)
        for idx, cluster in enumerate(public_clusters):
            public_clusters[idx] = np.delete(cluster, [ i for i,v in enumerate(cluster) if repertoire.clones.iloc[v]['frequency'] < minCloneCount] )
        public_clusters = [i for i in public_clusters if len(i) and len(i)>=minNodeCount ]
        cluster_candidates = [e for i in public_clusters for e in i]
        filteredClonotypes = prepared_clones.iloc[cluster_candidates] 
        df_list.append(filteredClonotypes)
        dictLen = len(dictIdx)
        for i,j in enumerate(cluster_candidates):
            dictIdx[i+dictLen] = j
    
    # print(cluster_candidates)
    new_clonotypes = pd.concat(df_list)
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
    public_clusters = clusters[-20:]
    for idx, cluster in enumerate(public_clusters):
        public_clusters[idx] = np.delete(cluster, [ i for i,v in enumerate(cluster) if new_clonotypes.iloc[v]['frequency'] < minCloneCount] )
    public_clusters = [i for i in public_clusters if len(i) and len(i)>=minNodeCount ]
    if absoulute_public == True :
        num_of_samples = np.unique(new_clonotypes["sampleID"].to_numpy()).shape[0]
        public_clusters = [i for i in public_clusters if np.unique(new_clonotypes.iloc[i]["sampleID"].to_numpy()).shape[0] == num_of_samples]
    else:
        public_clusters = [i for i in public_clusters if np.unique(new_clonotypes.iloc[i]["sampleID"].to_numpy()).shape[0] >= minimum_coverage]
    for idx,cluster in enumerate(public_clusters):
            public_clusters[idx] = np.array([dictIdx[i] for i in cluster])

    return public_clusters

if __name__ == "__main__":
   print(skeleton_similarity(test_csv_strategy().input(path),absoulute_public=True))
