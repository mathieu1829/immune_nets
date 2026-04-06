import numpy as np
from immune_nets.entities import ImmuneRepertoire
from immune_nets.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from immune_nets.creation.distance.levenshtein import levenshteinDistance
import igraph as ig

from immune_nets.creation.utils.pathManager import pathManager
from immune_nets.factories import ImmuneRepertoireFactory

path = pathManager().testDataPath / "healthy_test_clonotypes_0.csv"


def publicSimilarity(repertoire, frequencyCutoff = None, top_k = 20, absoulutePublic=False, minCoverage = 2, minClusterSize=2):
    """
    Identify public T-cell clusters across all samples in a repertoire.

    This function constructs a network from all clones in the repertoire using
    Levenshtein distance between TCR sequences, then performs clustering with
    the Fast Greedy algorithm. The top-k largest clusters are retained and
    filtered based on minimum cluster size and sample coverage.

    Note:
        - Clustering is performed on the global network including all samples.
        - Clusters containing clones from fewer than `minCoverage` samples are excluded.
        - Frequency cutoff can be applied to ignore low-frequency clones
          (currently requires assignment to take effect).

    Parameters
    ----------
    repertoire : ImmuneRepertoire
        The immune repertoire object containing clones with TCR sequences.
    frequencyCutoff : float, optional
        Minimum clone frequency to include in analysis. Clones below this frequency are excluded.
        Default is None (no cutoff).
    top_k : int, optional
        Number of largest clusters to retain. Default is 20.
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
    """

    prepared_clones = repertoire.clones.dropna(subset = ["tcra_aa", "tcrb_aa"]) 
    prepared_clones.name = repertoire.clones.name

    # If specified take into account only clones above frequency cutoff
    if frequencyCutoff is not None:
        prepared_clones[prepared_clones["frequency"] > frequencyCutoff]

    preparedRepertoire = ImmuneRepertoire(name="test repertoire" ,clones=prepared_clones)
    # print(prepared_clones)

    
    immuneNet = simpleBetaDistance(repertoire = preparedRepertoire,distance = levenshteinDistance(group = True),threshold = 2)
    df_net = immuneNet.graph

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
