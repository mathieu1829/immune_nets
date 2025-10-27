import optuna
from pathlib import Path
from itertools import combinations

from src.analysis.networkDistanceOptimization import objectiveBuilder

from src.models import Repertoire,Dataset
from src.db import engine
from sqlalchemy.orm import Session, selectinload
from sqlalchemy import insert,select,delete
from src.analysis.statDistance import EuclideanStatDistance
from src.mappers import RepertoireMapper

import time

groups = ["leukemia", "covid", "healthy"]
dataset_repertoires = {}
with Session(engine) as session: 
    for group in groups :
        stmt = select(Dataset).options(selectinload(Dataset.repertoires).selectinload(Repertoire.clonotypes)).where(Dataset.name == f"{group} test dataset")
        result = session.execute(stmt)
        dataset = result.scalars().first()
        dataset_repertoires[group] = [ RepertoireMapper.toImmuneRepertoire(repertoire) for repertoire in dataset.repertoires]

all_repertoires = { f"healthy vs {group}":{"healthy":dataset_repertoires["healthy"], group:dataset_repertoires[group]} for group in groups if not group == "healthy"}
all_repertoires["universal"] = dataset_repertoires

if __name__ == '__main__':
    for repertoire_group in all_repertoires:
        print(f"Running study for {repertoire_group} repertoires")
        analyzed_repertoires = all_repertoires[repertoire_group]
        study = optuna.create_study(direction="maximize")
        study.optimize(objectiveBuilder(analyzed_repertoires, EuclideanStatDistance()), n_trials=3)

        # Best result
        print("Best score:", study.best_value)
        print("Best params:", study.best_params)

