import optuna
from src.creation.io_strategies.test_csv_strategy import *
from pathlib import Path
from itertools import combinations

from src.analysis.networkDistanceOptimization import objectiveBuilder

# TO DO - get groups from db
groups = ["leukemia", "covid", "healthy"]

root_dir = Path(__file__).parent.parent.parent

leukemia_path = root_dir  / "tests/test_data/leukemia_test_clonotypes.csv" # leukemia
covid_path = root_dir / "tests/test_data/covid_test_clonotypes.csv" # covid
healthy_path = root_dir / "tests/test_data/healthy_test_clonotypes_1.csv" #healthy

# TO DO - get repertoires for each group from database
repertoire_list = [
        test_csv_strategy().input(leukemia_path),
        test_csv_strategy().input(covid_path),
        test_csv_strategy().input(healthy_path),
        ]
universal_repertoires = { group:[repertoire_list[i]]  for i,group in enumerate(groups)}

#gathering disease vs healthy datasets
all_repertoires = { f"healthy vs {group}":{"healthy":[repertoire_list[2]], group:[repertoire_list[i]]} for i, group in enumerate(groups) if not group == "healthy"}
#adding universal repertoires
all_repertoires["universal"] = universal_repertoires

if __name__ == '__main__':
    for repertoire_group in all_repertoires:
        print(f"Running study for {repertoire_group} repertoires")
        analyzed_repertoires = all_repertoires[repertoire_group]
        study = optuna.create_study(direction="maximize")
        study.optimize(objectiveBuilder(analyzed_repertoires), n_trials=20)

        # Best result
        print("Best score:", study.best_value)
        print("Best params:", study.best_params)


