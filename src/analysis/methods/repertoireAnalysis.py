import numpy as np
import pandas as pd
from src.creation.io_strategies.csv_strategy import csv_strategy
from src.creation.distance.hamming import hammingDistance 
from src.creation.algorithms.simple_beta_distance import simple_beta_distance
from src.creation.io_strategies.df_strategy import df_strategy
from src.creation.distance.negativeHamming import negativeHammingDistance
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.algorithms.common_methods import *
import igraph as ig
import uuid


df = pd.read_csv("10k_PBMC_5pv2_nextgem_Chromium_Controller_10k_PBMC_5pv2_nextgem_Chromium_Controller_vdj_t_clonotypes.csv").head(1000)
df.name = "aaa"

def repertoireAnalysis(repertoire, minNodeCount=10, minCloneCount=100):
    repertoire.clones['tcra_aa'] = repertoire.clones['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
    repertoire.clones['tcrb_aa'] = repertoire.clones['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
    tcra = np.unique(repertoire.clones['tcra_aa'].dropna().to_numpy())
    tcrb = np.unique(repertoire.clones['tcrb_aa'].dropna().to_numpy())
    all_tcr = np.unique(np.concatenate([tcra,tcrb]))

    num_of_tcra = tcra.shape[0]
    num_of_tcrb = tcrb.shape[0]
    num_of_all_tcr = all_tcr.shape[0]

    unique_tcra_distribution = np.array([ (tcra == seq).sum()  for seq in tcra ])
    unique_tcrb_distribution = np.array([ (tcrb == seq).sum()  for seq in tcrb ])
    unique_all_tcr_distribution = np.array([ (all_tcr == seq).sum()  for seq in all_tcr ])

    simpson_index_tcra = (unique_tcra_distribution**2).sum() 
    simpson_index_tcrb = (unique_tcrb_distribution**2).sum() 
    simpson_index_all_tcr = (unique_all_tcr_distribution**2).sum() 
    transform = lambda x : x * np.log(x)
    shannon_index_tcra = transform(unique_tcra_distribution).sum()
    shannon_index_tcrb = transform(unique_tcrb_distribution).sum()
    shannon_index_all_tcr = transform(unique_all_tcr_distribution).sum()

    ### if no graph as been provided create one
    ### public cluster creation
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
    public_count = sum([len(i) for i in public_clusters])
    public_cluster_count = len(public_clusters)


    immuneNet = simple_beta_distance(repertoire = repertoire,distance = negativeHammingDistance(),threshold = 1)
    df_net = immuneNet.network
    vertices = np.unique(df_net.to_numpy().flatten())
    net = ig.Graph(df_net.to_numpy())
    net.add_vertices(immuneNet.sampleSize - net.vcount())
    clusters = net.community_fastgreedy()
    clusters = list(clusters.as_clustering(clusters.optimal_count))
    clusters.sort(key = lambda x : len(x))
    private_clusters = clusters[-20:]
    for idx, cluster in enumerate(private_clusters):
        private_clusters[idx] = np.delete(cluster, [ i for i,v in enumerate(cluster) if repertoire.clones.iloc[v]['frequency'] <= minCloneCount] )
    private_clusters = [i for i in private_clusters if len(i) and len(i)<=minNodeCount ]
    private_count = sum([len(i) for i in private_clusters])
    private_cluster_count = len(private_clusters)
    
    
    return [num_of_tcra,num_of_tcrb,num_of_all_tcr, unique_tcra_distribution, unique_tcrb_distribution, unique_all_tcr_distribution, simpson_index_tcra, simpson_index_tcrb, simpson_index_all_tcr, shannon_index_tcra, shannon_index_tcrb, shannon_index_all_tcr, public_count, public_cluster_count, private_count, private_cluster_count]

if __name__ == "__main__":
    print(repertoireAnalysis(immuneRepertoire(df, {uuid.uuid4().hex:len(df)})))
