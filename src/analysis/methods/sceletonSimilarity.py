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
    public_clusters = clusters[-20:]
    for idx, cluster in enumerate(public_clusters):
        public_clusters[idx] = np.delete(cluster, [ i for i,v in enumerate(cluster) if repertoire.clones.iloc[v]['frequency'] >= minCloneCount] )
    public_clusters = [i for i in public_clusters if len(i) and len(i)>=minNodeCount ]
    
    cluster_candidates = [e for i in public_clusters for e in i]
    new_repertoire = immuneRepertoire(clones=repertoire.clones.iloc[cluster_candidates],sampleIDs=repertoire.sampleIDs) 
    immuneNet = simple_beta_distance(repertoire = new_repertoire,distance = hammingDistance(),threshold = -2)
    kdf_net = immuneNet.network
    print(df_net)
    vertices = np.unique(df_net.to_numpy().flatten())
    net = ig.Graph(df_net.to_numpy())
    net.add_vertices(immuneNet.sampleSize - net.vcount())
    clusters = net.community_fastgreedy()
    clusters = list(clusters.as_clustering(clusters.optimal_count))
    clusters.sort(key = lambda x : len(x))
    public_clusters = clusters[-20:]
    for idx, cluster in enumerate(public_clusters):
        public_clusters[idx] = np.delete(cluster, [ i for i,v in enumerate(cluster) if repertoire.clones.iloc[v]['frequency'] >= minCloneCount] )
    public_clusters = [i for i in public_clusters if len(i) and len(i)>=minNodeCount ]
    return public_clusters