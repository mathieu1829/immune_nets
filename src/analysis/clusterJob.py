import optuna
from pathlib import Path
from itertools import combinations

from src.analysis.networkDistanceOptimization import objectiveBuilder

from src.orm.models import Repertoire,Dataset
from src.orm.db import engine
from sqlalchemy.orm import Session, selectinload
from sqlalchemy import insert,select,delete

# TO DO - get groups from db
# groups = ["leukemia", "covid", "healthy"]
#
# root_dir = Path(__file__).parent.parent.parent
#
# leukemia_path = root_dir  / "tests/test_data/leukemia_test_clonotypes.csv" # leukemia
# covid_path = root_dir / "tests/test_data/covid_test_clonotypes.csv" # covid
# healthy_path = root_dir / "tests/test_data/healthy_test_clonotypes_1.csv" #healthy
#
# # TO DO - get repertoires for each group from database
# repertoire_list = [
#         Repertoire.from_csv(name="leukemia",desc=" ",path=leukemia_path),
#         Repertoire.from_csv(name="covid",desc=" ",path=covid_path),
#         Repertoire.from_csv(name="healthy",desc=" ",path=healthy_path),
#         ]
# universal_repertoires = { group:[repertoire_list[i]]  for i,group in enumerate(groups)}
#
# #gathering disease vs healthy datasets
# all_repertoires = { f"healthy vs {group}":{"healthy":[repertoire_list[2]], group:[repertoire_list[i]]} for i, group in enumerate(groups) if not group == "healthy"}
# #adding universal repertoires
# all_repertoires["universal"] = universal_repertoires
# all_repertoires = {}
#
with Session(engine) as session: 
    stmt1 = select(Dataset).options(selectinload(Dataset.repertoires).selectinload(Repertoire.clonotypes)).where(Dataset.name == "healthy")
    result1 = session.execute(stmt1)
    # print(list(result1.all()))
    healthy_dataset: Dataset | None = result1.scalars().first()

    stmt2 = select(Dataset).options(selectinload(Dataset.repertoires).selectinload(Repertoire.clonotypes)).where(Dataset.name == "leukemia")
    result2 = session.execute(stmt2)
    # print(list(result2.all()))
    leukemia_dataset = result2.scalars().first()

    if leukemia_dataset is None or healthy_dataset is None:
        raise ValueError("One of the repertoires was not found")
    all_repertoires = {"healthy vs leukemia": {"healthy":healthy_dataset.repertoires , "leukemia":leukemia_dataset.repertoires}}

if __name__ == '__main__':
    for repertoire_group in all_repertoires:
        if repertoire_group != "healthy vs leukemia":
            continue
        print(f"Running study for {repertoire_group} repertoires")
        analyzed_repertoires = all_repertoires[repertoire_group]
        study = optuna.create_study(direction="maximize")
        study.optimize(objectiveBuilder(analyzed_repertoires), n_trials=3)

        # Best result
        print("Best score:", study.best_value)
        print("Best params:", study.best_params)


