import optuna
import numpy as np
from itertools import combinations
from multiprocessing import Pool

from .computeGraphStats import computeGraphStats 

from src.creation.distance.alignment import sequenceAligner
from src.creation.distance.levenshtein import levenshteinDistance 

from src.analysis.scoringParadigms import ScoringParadigm
from src.cluster.utils.commonMethods import createTestGroups

def optimizerObjectiveBuilder(repertoireDatasets, scoringParadigm: ScoringParadigm, rank):
    scoringParadigmFun = scoringParadigm.compute_score
    def objective(trial):
        threshold = trial.suggest_float("threshold",low=0.2,high=0.4)
        distance = trial.suggest_categorical("distance", ["alignment", "levenshtein"])
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
            case _:
                distance_fun = levenshteinDistance()

        repertoireStats = { group:[None for _ in repertoireDatasets[group]] for group in repertoireDatasets}
        args = [{"repertoire": repertoire,
                 "repertoireDataset": repertoireDataset,
                 "repertoireIdx": idx,
                 "algorithm_name": algorithm_name,
                 "threshold": threshold,
                 "distance_fun": distance_fun,
                 "worldRank": rank,
                 } for repertoireDataset in repertoireDatasets for idx, repertoire in enumerate(repertoireDatasets[repertoireDataset]) ]

        with Pool(3) as pool:
            for result in pool.imap_unordered(computeGraphStats, args, chunksize=1):
                if result is None:
                    pool.terminate()
                    return 0.0
                else:
                    repertoireStats[result["dataset"]][result["idx"]] = result["value"]

        testGroups = createTestGroups(repertoireStats)
        args = [testGroup for testGroup in testGroups.values()]
        with Pool(3) as p:
            resultList = p.map(scoringParadigmFun, args)

        between = resultList[0]
        within = np.mean(resultList[1:])
        trial.set_user_attr(f"between", between)
        trial.set_user_attr(f"within", within)
        if between == 0.0:
            return 0.0

        # within * 0.95
        # within * heuristic progress (temperature)
        return ((between - within)/between)  
    return objective









