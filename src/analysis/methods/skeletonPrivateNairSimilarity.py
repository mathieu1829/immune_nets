import numpy as np
import pandas as pd
from src.creation.immuneRepertoire import ImmuneRepertoire
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.distance.levenshtein import levenshteinDistance
from src.analysis.methods.skeletonPublicNairSimilarity import skeletonPublicNairSimilarity
import igraph as ig
import itertools

from src.creation.utils.pathManager import pathManager
from src.factories import ImmuneRepertoireFactory

path = pathManager().testDataPath / "healthy_test_clonotypes_0.csv"


def skeletonPrivateNairSimilarity(repertoire, top_k = 20, absoulutePublic=False, minCoverage = 2, minClusterSize=2, minPrivateClusterSize=1):
    public_clusters = skeletonPublicNairSimilarity(repertoire=repertoire, 
                                       top_k=top_k,absoulutePublic=absoulutePublic,
                                       minCoverage=minCoverage,
                                       minClusterSize=minClusterSize)
    public_clusters = set(list(itertools.chain(*public_clusters)))

    prepared_clones = repertoire.clones.dropna(subset = ["tcra_aa", "tcrb_aa"]) 
    # print(prepared_clones)

    globalPrivateClusters = []
    # Creating network and cluster analysis for each sample
    for sampleID in np.unique(prepared_clones["sampleID"]):
        # print(f"sampleID: {sampleID}")
        #creating network for selected sample
        sampleClones = prepared_clones.iloc[ prepared_clones["sampleID"].to_numpy() == sampleID]
        sampleClones.name = repertoire.clones.name
        sampleRepertoire = ImmuneRepertoire(name="test repertoire", clones=sampleClones)
        immuneNet = simpleBetaDistance(repertoire = sampleRepertoire,distance = levenshteinDistance(group = True),threshold = 2)
        df_net = immuneNet.graph

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
        private_clusters = clusters[-top_k:]
        # print(private_clusters)

        # assigning indices from original df
        private_clusters = [ sampleClones.iloc[cluster].index.tolist() for cluster in private_clusters]

        # collecting clusters that have no overlap with public and have big enough size
        samplePrivateClusters = [ cluster for cluster in private_clusters if not set(cluster) & public_clusters and len(cluster) >= minPrivateClusterSize]
        globalPrivateClusters.append(samplePrivateClusters)

    return globalPrivateClusters 


   

if __name__ == "__main__":
   print(skeletonPrivateNairSimilarity(ImmuneRepertoireFactory.fromCSVTest(path)))
