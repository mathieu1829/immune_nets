import optuna
import numpy as np
from itertools import combinations

from src.entities import GraphStats
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance

def compareGroups(repertoires, algorithm_name, threshold, distance_fun, scoringParadigmFun):
    match algorithm_name:
        case "simpleBetaDistance":
            algorithm = simpleBetaDistance
        case "simpleVectorBetaDistance":
            algorithm = simpleVectorBetaDistance
        case _:
            algorithm = simpleBetaDistance #default

    groups = [group for group in repertoires]
    group_results = { group:[] for group in groups}
    for group in groups:
        for repertoire in repertoires[group]:
            network = algorithm(
                        repertoire=repertoire,
                        distance=distance_fun,
                        threshold=threshold
                    )
            stats = GraphStats(network)
            if stats.numEdges == 0:
                return 0.0
            V = stats.verticeNum
            if stats.numEdges == ((V*(V-1))/2):
                return 0.0
            group_results[group].append(stats)
            # print(f"{group} stats: {str(stats.toList())}")
                
    return scoringParadigmFun(group_results) 

