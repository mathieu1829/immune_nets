import optuna
import numpy as np
from itertools import combinations
from multiprocessing import Pool


from immune_nets.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from immune_nets.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from immune_nets.creation.distance.alignment import sequenceAligner
from immune_nets.creation.distance.levenshtein import levenshteinDistance 
from immune_nets.analysis.scoringParadigms import ScoringParadigm
from .utils.computeGraphStats import computeGraphStats

class initialAnalysisObjectiveBuilder:
    """
    
    """
    def __init__(self,repertoires, scoringParadim: ScoringParadigm, rank: int, statComputingPoolSize=1):
        self.repertoires = repertoires
        self.scoringParadigmFun = scoringParadim.compute_score
        self.rank = rank
        self.statComputingPoolSize = statComputingPoolSize

    def __call__(self, trial):
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

        repertoireStats = { group:[None for _ in self.repertoires[group]] for group in self.repertoires}
        args = [{"repertoire": repertoire,
                 "repertoireDatasetName": testGroup,
                 "repertoireIdx": idx,
                 "algorithm_name": algorithm_name,
                 "threshold": threshold,
                 "distance_fun": distance_fun,
                 "worldRank": self.rank,
                 } for testGroup in self.repertoires for idx, repertoire in enumerate(self.repertoires[testGroup]) ]

        with Pool(self.statComputingPoolSize) as pool:
            for result in pool.imap_unordered(computeGraphStats, args, chunksize=1):
                if result is None:
                    pool.terminate()
                    return 0.0
                else:
                    repertoireStats[result["dataset"]][result["idx"]] = result["value"]

        result = self.scoringParadigmFun(repertoireStats)
                    
        return result









