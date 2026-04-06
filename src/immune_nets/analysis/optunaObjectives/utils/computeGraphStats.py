import optuna
import numpy as np
from itertools import combinations

from immune_nets.entities import GraphStats
from immune_nets.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from immune_nets.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance

def computeGraphStats(args):
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

