import optuna
import numpy as np
from itertools import combinations
from mpi4py import MPI

from .compareGroups import compareGroups

from src.creation.distance.alignment import sequenceAligner
from src.creation.distance.levenshtein import levenshteinDistance 
from src.analysis.scoringParadigms import ScoringParadigm

def optimizerObjectiveBuilder(repertoires, repertoire_group, scoringParadigm: ScoringParadigm, cluster: MPI.Comm):
    scoringParadigmFun = scoringParadigm.compute_score
    def objective(trial):
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



        for i in range(1,cluster.Get_size()):
            cluster.send(obj=True,dest=i, tag=1)
            cluster.send(obj=repertoires, dest=i, tag=1) 
            cluster.send(obj=algorithm_name, dest=i, tag=1) 
            cluster.send(obj=threshold, dest=i, tag=1) 
            cluster.send(obj=distance_fun, dest=i, tag=1) 
            cluster.send(obj=scoringParadigmFun, dest=i, tag=1) 

        result = compareGroups(repertoires,
                               algorithm_name,
                               threshold,
                               distance_fun,
                               scoringParadigmFun) 

        groupList = cluster.gather(repertoire_group, root=0)
        resultList = cluster.gather(result, root=0)

        between = 0.0
        within = []
        for result, repertoire_group_name in zip(resultList, groupList):
            if not "_" in repertoire_group_name:
                between = result
            else:
                within.append(result)

        within = np.mean(within)
        if between == 0.0:
            return 0.0

        # within * 0.95
        # within * heuristic progress (temperature)
        return ((between - 0.8 * within)/between)  
    return objective









