import optuna
import numpy as np
from itertools import combinations

from src.analysis.methods.graphStats import GraphStats
from src.mappers import GraphStatsMapper
from src.creation.algorithms.simpleDistance import simpleDistance
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from src.creation.distance.alignment import sequenceAligner
from src.creation.distance.levenshtein import levenshteinDistance 
from scipy.spatial.distance import euclidean

def objectiveBuilder(repertoires):
    groups = [group for group in repertoires]
    def objective(trial):
        group_results = { group:[] for group in groups}
        threshold = trial.suggest_float("threshold",low=0.2,high=0.4)
        distance = trial.suggest_categorical("distance", ["alignment", "levenshtein"])
        distance_fun = None
        algorithm_name = trial.suggest_categorical("algorithm_name", ["simpleBetaDistance", "simpleVectorBetaDistance"])
        algorithm = None
        match distance:
            case "alignment":
                substitution_matrix = trial.suggest_categorical("substitution_matrix", [
                    "PAM250",
                    "PAM30",
                    "PAM70",
                    "BLOSUM45",
                    "BLOSUM50",
                    "BLOSUM62",
                    "BLOSUM80",
                    "BLOSUM90"
                    ])
                distance_fun = sequenceAligner(substitution_matrix) 
            case "levenshtein":
                distance_fun = levenshteinDistance()

        match algorithm_name:
            case "simpleBetaDistance":
                algorithm = simpleBetaDistance
            case "simpleVectorBetaDistance":
                algorithm = simpleVectorBetaDistance


        for group in groups:
            for repertoire in repertoires[group]:
                network = algorithm(
                            repertoire=repertoire,
                            distance=distance_fun,
                            threshold=threshold
                        )
                stats = GraphStats(network)
                group_results[group].append(GraphStatsMapper.toList(stats))
                # print(f"{group} stats: {str(stats.toList())}")
                    
        inter_group_distances = []
        for combo in combinations(groups,2):
            group_a = group_results[combo[0]]
            group_b = group_results[combo[1]]
            inter_group_distance = np.array([euclidean(a,b) for a in group_a for b in group_b ])
            inter_group_distances.append(inter_group_distance.mean())
            # print(f"distance between group {combo[0]} and {combo[1]} is {inter_group_distance.mean()}")
        
        inter_group_distances = np.array(inter_group_distances)
        # print(f"mean: {inter_group_distances.mean()}")
        # print(f"std: {np.std(inter_group_distances)}")
        # print(f"var: {np.var(inter_group_distances)}")
        return inter_group_distances.mean() - np.std(inter_group_distances) 
    return objective









