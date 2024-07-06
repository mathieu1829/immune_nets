import numpy as np
import pandas as pd
from src.creation.io_strategies.csv_strategy import csv_strategy
from src.creation.distance.hamming import hammingDistance 
from src.creation.distance.groupHamming import groupHammingDistance 
from src.creation.distance.groupLevenshtein import groupLevenshteinDistance 
from src.creation.algorithms.simple_beta_distance import simple_beta_distance
from src.creation.io_strategies.df_strategy import df_strategy
from src.creation.distance.negativeHamming import negativeHammingDistance
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.algorithms.common_methods import *
import igraph as ig
import uuid
from src.creation.utils.pathManager import pathManager
from src.creation.io_strategies.test_csv_strategy import *

path = pathManager().testDataPath / "bigTest.csv"

class repertoireAnalysis:
    def __init__(self,repertoire, minNodeCount=10, minCloneCount=100):
        repertoire.clones['tcra_aa'] = repertoire.clones['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        repertoire.clones['tcrb_aa'] = repertoire.clones['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        tcra = np.unique(repertoire.clones['tcra_aa'].dropna().to_numpy())
        tcrb = np.unique(repertoire.clones['tcrb_aa'].dropna().to_numpy())
        all_tcr = np.unique(np.concatenate([tcra,tcrb]))

        self.num_of_tcra = tcra.shape[0]
        self.num_of_tcrb = tcrb.shape[0]
        self.num_of_all_tcr = all_tcr.shape[0]

        self.unique_tcra_distribution = np.array([ (tcra == seq).sum()  for seq in tcra ])
        self.unique_tcrb_distribution = np.array([ (tcrb == seq).sum()  for seq in tcrb ])
        self.unique_all_tcr_distribution = np.array([ (all_tcr == seq).sum()  for seq in all_tcr ])

        self.simpson_index_tcra = (self.unique_tcra_distribution**2).sum() 
        self.simpson_index_tcrb = (self.unique_tcrb_distribution**2).sum() 
        self.simpson_index_all_tcr = (self.unique_all_tcr_distribution**2).sum() 
        transform = lambda x : x * np.log(x)
        self.shannon_index_tcra = transform(self.unique_tcra_distribution).sum()
        self.shannon_index_tcrb = transform(self.unique_tcrb_distribution).sum()
        self.shannon_index_all_tcr = transform(self.unique_all_tcr_distribution).sum()

        ### if no graph as been provided create one
        ### public cluster creation
        immuneNet = simple_beta_distance(repertoire = repertoire,distance = groupLevenshteinDistance(),threshold = 2)
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
        self.public_count = sum([len(i) for i in public_clusters])
        self.public_cluster_count = len(public_clusters)


        
    def toList(self): 
        return [self.num_of_tcra, self.num_of_tcrb, self.num_of_all_tcr, self.unique_tcra_distribution, self.unique_tcrb_distribution, self.unique_all_tcr_distribution, self.simpson_index_tcra, self.simpson_index_tcrb, self.simpson_index_all_tcr, self.shannon_index_tcra, self.shannon_index_tcrb, self.shannon_index_all_tcr, self.public_count, self.public_cluster_count]

if __name__ == "__main__":
    print(repertoireAnalysis(test_csv_strategy().input(path)))
