from optuna import Trial
import numpy as np
from multiprocessing import Pool

from .utils.computeGraphStats import computeGraphStats 
from .objectiveBuilder import objectiveBuilder

from immune_nets.creation.distance.alignment import sequenceAligner
from immune_nets.creation.distance.levenshtein import levenshteinDistance 

from immune_nets.analysis.scoringParadigms import ScoringParadigm
from immune_nets.cluster.utils.commonMethods import createTestGroups
from immune_nets.entities.immuneRepertoire import ImmuneRepertoire

class OptimizerObjectiveBuilder(objectiveBuilder):
    """
    Objective builder for optimizing separation between groups of similarity networks.

    Constructs an objective function that:
    - generates similarity networks for all repertoires,
    - partitions them into test groups and cohorts,
    - maximizes distances between cohorts from different disease states,
    - minimizes distances within cohorts from the same disease state.

    The objective encourages parameter configurations that produce well-separated
    disease groups while maintaining internal consistency within each group. In effect
    it is trying to find network construction parameters which give maximal distance between 
    groups and minimal distance within groups.

    If any generated network is empty or fully connected, the evaluation is aborted
    and a value of 0.0 is returned.
    
    The final objective value is computed as:

    (between - within) / between

    where:
    - between: distance between cohorts of different disease states,
    - within: mean distance within cohorts of the same disease state.
    """
    def __init__(self, repertoireDatasets: dict[str, list[ImmuneRepertoire]], scoringParadigm: ScoringParadigm, rank: int, statComputingPoolSize: int, resultGatheringPoolSize: int):
        """
        :param repertoireDatasets: Dictionary mapping disease state names to lists of ``ImmuneRepertoire`` objects.
        :param scoringParadigm: Scoring paradigm used to find the distance between cohorts of networks
        :param rank: Identifier of the process (used for parallel execution and logging).
        :param statComputingPoolSize: number of processes which concurently compute similarity networks from repertoires.
        :param resultGatheringPoolSize: number of processes which compute distances of individual test groups that will be later used to compute final value of returned by objective function.
        """

        self.repertoireDatasets = repertoireDatasets
        self.scoringParadigmFun = scoringParadigm.compute_score
        self.rank = rank
        if statComputingPoolSize > 0:
            self.statComputingPoolSize = statComputingPoolSize
        else:
            raise ValueError("Number of threads in stat computing pool (statComputingPoolSize) must be higher than 0")

        if resultGatheringPoolSize > 0:
            self.resultGatheringPoolSize = resultGatheringPoolSize
        else:
            raise ValueError("Number of threads in result gathering pool (resultGatheringPoolSize) must be higher than 0")


    def __call__(self, trial: Trial) -> float:
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

        repertoireStats = { group:[None for _ in self.repertoireDatasets[group]] for group in self.repertoireDatasets}
        args = [{"repertoire": repertoire,
                 "repertoireDatasetName": repertoireDatasetName,
                 "repertoireIdx": idx,
                 "algorithm_name": algorithm_name,
                 "threshold": threshold,
                 "distance_fun": distance_fun,
                 "worldRank": self.rank,
                 } for repertoireDatasetName in self.repertoireDatasets for idx, repertoire in enumerate(self.repertoireDatasets[repertoireDatasetName]) ]

        with Pool(self.statComputingPoolSize) as pool:
            for result in pool.imap_unordered(computeGraphStats, args, chunksize=1):
                if result is None:
                    pool.terminate()
                    return 0.0
                else:
                    repertoireStats[result["dataset"]][result["idx"]] = result["value"]

        testGroups = createTestGroups(repertoireStats)
        args = [testGroup for testGroup in testGroups.values()]
        with Pool(self.resultGatheringPoolSize) as p:
            resultList = p.map(self.scoringParadigmFun, args)

        between = resultList[0]
        within = np.mean(resultList[1:])
        trial.set_user_attr(f"between", between)
        trial.set_user_attr(f"within", within)
        if between == 0.0:
            return 0.0

        # within * 0.95
        # within * heuristic progress (temperature)
        return float((between - within)/between)  









