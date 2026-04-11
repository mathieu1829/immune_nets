import optuna
import numpy as np
from itertools import combinations

from immune_nets.entities import GraphStats
from immune_nets.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from immune_nets.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance

def computeGraphStats(args):
    """
    Compute a similarity network and its statistics for a given repertoire.

    This function is designed for use with ``multiprocessing.Pool``, where
    arguments must be passed as a single object.

    :param args: Dictionary containing the following keys:
        - ``repertoire`` (ImmuneRepertoire): Input repertoire.
        - ``repertoireDatasetName`` (str): Name of the dataset/cohort.
        - ``repertoireIdx`` (int): Index of the repertoire within the dataset.
        - ``algorithm_name`` (str): Name of the network construction algorithm.
        - ``threshold`` (float): Threshold parameter for network construction.
        - ``distance_fun`` (callable): Distance function used by the algorithm.
        - ``worldRank`` (int): Identifier of the parent process (for logging).
    :return: Dictionary with computed statistics (``value``, ``dataset``, ``idx``),
        or ``None`` if the generated network is empty or fully connected.
    """
    repertoire = args["repertoire"]
    repertoireDatasetName = args["repertoireDatasetName"]
    repertoireIdx = args["repertoireIdx"]
    algorithm_name = args["algorithm_name"]
    threshold = args["threshold"]
    distance_fun = args["distance_fun"]
    worldRank = args["worldRank"]
    
    match algorithm_name:
        case "simpleBetaDistance":
            algorithm = simpleBetaDistance
        case "simpleVectorBetaDistance":
            algorithm = simpleVectorBetaDistance
        case _:
            algorithm = simpleBetaDistance #default

    print(f"Worker of the main process {worldRank} is commencing computation for repertoire {repertoireIdx} of group: {repertoireDatasetName}")


    network = algorithm(
                repertoire=repertoire,
                distance=distance_fun,
                threshold=threshold
            )
    stats = GraphStats(network)

    #Maybe change for within vs between
    if stats.numEdges == 0:
        return None
    V = stats.verticeNum
    if stats.numEdges == ((V*(V-1))/2):
        return None
    # print(f"{group} stats: {str(stats.toList())}")

    print(f"Worker of the main process {worldRank} has finished processing {repertoireDatasetName} ")
                
    return {"value":stats, "dataset":repertoireDatasetName, "idx": repertoireIdx} 

