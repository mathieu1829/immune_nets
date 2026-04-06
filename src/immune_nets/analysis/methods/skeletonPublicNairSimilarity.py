import numpy as np
import pandas as pd
from immune_nets.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from immune_nets.creation.distance.levenshtein import levenshteinDistance
import igraph as ig

from immune_nets.creation.utils.pathManager import pathManager
from immune_nets.entities import ImmuneRepertoire
from immune_nets.factories import ImmuneRepertoireFactory

path = pathManager().testDataPath / "healthy_test_clonotypes_0.csv"


def skeletonPublicNairSimilarity(repertoire, top_k = 20, absoulutePublic=False, minCoverage = 2, minClusterSize=2):
    """
    Identify public T-cell clusters using a skeleton-based per-sample approach as shown in NAIR artice.

    This function first constructs separate networks per sample, identifies
    the top-k clusters per sample, and selects a representative "skeleton" clone
    from each cluster (the highest-frequency clone). Skeleton clones from all
    samples are then merged into a global network, clustered, and filtered
    based on sample coverage and minimum cluster size. The final clusters
    are expanded to include all original clones from each representative's
    initial cluster.

    Key differences from `publicSimilarity`:
        - Uses a per-sample skeletonization step before global clustering.
        - Reduces cross-sample noise and computational load.
        - Expands skeleton clusters back to include all original clones.

    Parameters
    ----------
    repertoire : ImmuneRepertoire
        The immune repertoire object containing clones with TCR sequences.
    top_k : int, optional
        Number of largest clusters to retain per sample. Default is 20.
    absoulutePublic : bool, optional
        If True, considers clusters as public only if they are present in all samples.
        Default is False.
    minCoverage : int, optional
        Minimum number of samples in which a cluster must appear to be retained.
        Default is 2.
    minClusterSize : int, optional
        Minimum size (number of clones) of clusters to retain. Default is 2.

    Returns
    -------
    list[list[int]]
        List of public clusters. Each cluster is a list of indices referring
        to the clones in the original repertoire DataFrame.

    Notes
    -----
    - Skeletonization involves selecting a representative clone per cluster per sample.
    - Cluster expansion ensures that each final cluster contains all original clones
      from the representative clone's initial cluster.
    """

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
    skeleton_repertoire = ImmuneRepertoire(name="test repertoire", clones=skeleton_clones)
    immuneNet = simpleBetaDistance(repertoire = skeleton_repertoire,distance = levenshteinDistance(group = True),threshold = 2)
    df_net = immuneNet.graph

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
   print(skeletonPublicNairSimilarity(ImmuneRepertoireFactory.fromCSVTest(path),absoulutePublic=True))
