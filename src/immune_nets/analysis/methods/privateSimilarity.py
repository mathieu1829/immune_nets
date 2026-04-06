import numpy as np
import igraph as ig
import itertools

from immune_nets.entities import ImmuneRepertoire
from immune_nets.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from immune_nets.creation.distance.levenshtein import levenshteinDistance
from immune_nets.analysis.methods.publicSimilarity import publicSimilarity

from immune_nets.creation.utils.pathManager import pathManager
from immune_nets.factories import ImmuneRepertoireFactory

path = pathManager().testDataPath / "healthy_test_clonotypes_0.csv"


def privateSimilarity(repertoire, frequencyCutoff = None, top_k = 20, absoulutePublic=False, minCoverage = 2, minClusterSize=2, minPrivateClusterSize=1):
    """
    Identify private T-cell clusters in a repertoire by constructing a global network.

    This function first identifies public clusters using the `publicSimilarity` method,
    then builds a network across all samples in the repertoire based on
    Levenshtein distances between TCR sequences. Clusters are extracted from the network,
    filtered to remove overlaps with public clusters, and grouped by sample.

    Note:
        - Clustering is performed globally across all samples.
        - Frequency cutoff filter (if provided) is applied to clone frequencies.
        - Clusters containing clones from multiple samples are excluded.

    Parameters
    ----------
    repertoire : ImmuneRepertoire
        The immune repertoire object containing clones with TCR sequences.
    frequencyCutoff : float, optional
        Minimum clone frequency to include in analysis. Clones below this frequency are excluded.
        Default is None (no cutoff).
    top_k : int, optional
        Number of largest clusters to retain per sample. Default is 20.
    absoulutePublic : bool, optional
        Whether to use absolute public clonotype definition in `publicSimilarity`.
        Default is False.
    minCoverage : int, optional
        Minimum coverage for public clusters. Default is 2.
    minClusterSize : int, optional
        Minimum size of clusters to consider in public similarity. Default is 2.
    minPrivateClusterSize : int, optional
        Minimum size of private clusters to retain. Default is 1.

    Returns
    -------
    list[list[list[int]]]
        Nested list of private clusters grouped by sample.
        Outer list: samples
        Middle list: clusters in sample
        Inner list: indices of clones in the original DataFrame
    """
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

    preparedRepertoire = ImmuneRepertoire(name="test repertoire",clones=prepared_clones)
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
   print(privateSimilarity(ImmuneRepertoireFactory.fromCSVTest(path)))
