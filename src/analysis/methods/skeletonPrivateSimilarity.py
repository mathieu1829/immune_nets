import numpy as np
import pandas as pd
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.algorithms.simple_beta_distance import simple_beta_distance
from src.creation.distance.hamming import hammingDistance
import igraph as ig

def skeleton_similarity(repertoire, minNodeCount=10, minCloneCount=100):
    immuneNet = simple_beta_distance(repertoire = repertoire,distance = hammingDistance(),threshold = -2)
    df_net = immuneNet.network
    print(df_net)
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
    new_repertoire = immuneRepertoire(clones=repertoire.clones.iloc[cluster_candidates],sampleIDs=repertoire.sampleIDs[cluster_candidates]) 
    immuneNet = simple_beta_distance(repertoire = new_repertoire,distance = hammingDistance(),threshold = -2)
    df_net = immuneNet.network
    print(df_net)
    vertices = np.unique(df_net.to_numpy().flatten())
    net = ig.Graph(df_net.to_numpy())
    net.add_vertices(immuneNet.sampleSize - net.vcount())
    clusters = net.community_fastgreedy()
    clusters = list(clusters.as_clustering(clusters.optimal_count))
    clusters.sort(key = lambda x : len(x))
    private_clusters = clusters
    for idx, cluster in enumerate(private_clusters):
        private_clusters[idx] = np.delete(cluster, [ i for i,v in enumerate(cluster) if repertoire.clones.iloc[v]['frequency'] >= minCloneCount] )
    private_clusters = [i for i in private_clusters if len(i) and np.unique(new_repertoire.sampleIDs[i]) == 1 ]
    return private_clusters
