import optuna
import numpy as np
from itertools import combinations

from src.analysis.methods.graphlet_composition import graphletComposition
from src.creation.algorithms.simple_distance import simple_distance
from src.creation.distance.alignment import sequenceAligner
from src.creation.distance.levenshtein import levenshteinDistance 
from scipy.spatial.distance import euclidean
from src.creation.io_strategies.test_csv_strategy import *

# TO DO - get groups from db
groups = ["leukemia", "covid", "healthy"]


leukemia_path = "/home/myc0plasmus/Downloads/T_PLL_sorted_5pv2_nextgem_vdj_t_clonotypes.csv" # leukemia
covid_path = "/home/myc0plasmus/Documents/studia/immune-nets-datasets/clonotypes.csv" # covid
healthy_path = "/home/myc0plasmus/Downloads/10k_BMMNC_5pv2_nextgem_intron_10k_BMMNC_5pv2_nextgem_intron_vdj_t_clonotypes.csv" #healthy

# TO DO - get repertoires for each group from database
repertoire_list = [
        test_csv_strategy().input(leukemia_path),
        test_csv_strategy().input(covid_path),
        test_csv_strategy().input(healthy_path),
        ]
repertoires = { group:[repertoire_list[i]]  for i,group in enumerate(groups)}


def objective(trial):
    group_results = { group:[] for group in groups}
    threshold = trial.suggest_float("threshold",low=0.01,high=0.5)
    distance = trial.suggest_categorical("distance", ["alignment", "levenshtein"])
    distance_fun = None
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

    for group in groups:
        for repertoire in repertoires[group]:
            network = simple_distance(
                        repertoire=repertoire,
                        distance=distance_fun,
                        threshold=threshold
                    )
            try:
                stats = graphletComposition(network)
            except ZeroDivisionError: 
                return 0.0
            group_results[group].append(stats.toList())
                
    inter_group_distances = []
    for combo in combinations(groups,2):
        group_a = group_results[combo[0]]
        group_b = group_results[combo[1]]
        inter_group_distance = np.array([euclidean(a,b) for a in group_a for b in group_b ])
        inter_group_distances.append(inter_group_distance.mean())
    
    inter_group_distances = np.array(inter_group_distances)
    return inter_group_distances.mean() - np.var(inter_group_distances)

if __name__ == '__main__':
    study = optuna.create_study(direction="maximize")
    study.optimize(objective, n_trials=20)

    # Best result
    print("Best score:", study.best_value)
    print("Best params:", study.best_params)









