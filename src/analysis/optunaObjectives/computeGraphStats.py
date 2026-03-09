import optuna
import numpy as np
from itertools import combinations

from src.entities import GraphStats
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance

def computeGraphStats(args):
    repertoire = args["repertoire"]
    repertoireDataset = args["repertoireDataset"]
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

    print(f"Worker of the main process {worldRank} is commencing computation for repertoire {repertoireIdx} of group: {repertoireDataset}")


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

    print(f"Worker of the main process {worldRank} has finished processing {repertoireDataset} ")
                
    return {"value":stats, "dataset":repertoireDataset, "idx": repertoireIdx} 

