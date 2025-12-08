import optuna
import numpy as np
from itertools import combinations

from src.entities import GraphStats
from src.creation.algorithms.simpleDistance import simpleDistance
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from src.creation.distance.alignment import sequenceAligner
from src.creation.distance.levenshtein import levenshteinDistance 

from src.analysis.scoringParadigms import ScoringParadigm

def objectiveBuilder(repertoires, scoringParadim: ScoringParadigm):
    groups = [group for group in repertoires]
    scoringParadimFun = scoringParadim.compute_score
    def objective(trial):
        group_results = { group:[] for group in groups}
        threshold = trial.suggest_float("threshold",low=0.2,high=0.4)
        distance = trial.suggest_categorical("distance", ["alignment", "levenshtein"])
        distance_fun = None
        algorithm_name = trial.suggest_categorical("algorithm_name", ["simpleBetaDistance", "simpleVectorBetaDistance"])
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
            case _:
                algorithm = simpleBetaDistance #default


        for group in groups:
            for repertoire in repertoires[group]:
                network = algorithm(
                            repertoire=repertoire,
                            distance=distance_fun,
                            threshold=threshold
                        )
                stats = GraphStats(network)
                group_results[group].append(stats)
                # print(f"{group} stats: {str(stats.toList())}")
                    
        return scoringParadimFun(group_results) 
    return objective









