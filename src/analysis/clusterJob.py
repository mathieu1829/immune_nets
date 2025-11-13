import optuna
from pathlib import Path
from itertools import combinations
import numpy as np
import uuid
import pickle

from src.analysis.networkDistanceOptimization import objectiveBuilder

from src.models import Repertoire,Dataset
from src.db import engine
from sqlalchemy.orm import Session, selectinload
from sqlalchemy import insert,select,delete
from src.analysis.statDistance import WassersteinStatDistance
from src.mappers import RepertoireMapper, ImmuneNetworkMapper

from src.creation.immuneRepertoire import ImmuneRepertoire
from src.creation.immuneNetwork import ImmuneNetwork
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from src.creation.distance.levenshtein import levenshteinDistance
from src.creation.distance.alignment import sequenceAligner


with open("a.csv","w") as f:
    f.write("a")

def shuffleGroupAndDivide(groupedRepertoires, group):
  repertoireGroup = groupedRepertoires[group]
  np.random.shuffle(repertoireGroup)
  mid = len(repertoireGroup) // 2
  return {f"{group}_1":repertoireGroup[:mid], f"{group}_2":repertoireGroup[mid:]}

def makeBestNetwork(repertoire: ImmuneRepertoire, study) -> ImmuneNetwork:
    threshold = study.best_params["threshold"]
    distance = study.best_params["distance"]
    distance_fun = None
    algorithm_name = study.best_params["algorithm_name"]
    algorithm = None
    match distance:
        case "alignment":
            substitution_matrix = study.best_params["substitution_matrix"]
            distance_fun = sequenceAligner(substitution_matrix)
        case "levenshtein":
            distance_fun = levenshteinDistance()

    match algorithm_name:
        case "simpleBetaDistance":
            algorithm = simpleBetaDistance
        case "simpleVectorBetaDistance":
            algorithm = simpleVectorBetaDistance
        case _:
            algorithm = simpleBetaDistance


    network = algorithm(repertoire=repertoire, threshold=threshold, distance=distance_fun)
    return network
    

groups = ["covid","healthy"]
groupSize = [2,2]

dataset_repertoires = {}
with Session(engine) as session: 
    for group in groups :
        stmt = select(Dataset).options(selectinload(Dataset.repertoires).selectinload(Repertoire.clonotypes)).where(Dataset.name == f"{group} test dataset")
        result = session.execute(stmt)
        dataset = result.scalars().first()
        dataset_repertoires[group] = [ RepertoireMapper.toImmuneRepertoire(repertoire) for repertoire in dataset.repertoires]

all_repertoires = { f"healthy vs {group}":{"healthy":dataset_repertoires["healthy"], group:dataset_repertoires[group]} for group in groups if not group == "healthy"}
# all_repertoires["universal"] = dataset_repertoires
for group in groups:
  all_repertoires[f"{group}_1 vs {group}_2"] = shuffleGroupAndDivide(dataset_repertoires, group)

if __name__ == '__main__':
    distributionNames = ["degreeDistribution", "componentSizeDistribution", "componentProportionDistribution"]
    distributionCombos = list(combinations(distributionNames, 2))
    runId = uuid.uuid4()

    for combo in distributionCombos:
        results = {}
        for repertoire_group in all_repertoires:
            print(f"Running study for {repertoire_group} repertoires")
            analyzed_repertoires = all_repertoires[repertoire_group]
            study = optuna.create_study(direction="maximize")
            objectiveFunction = objectiveBuilder(repertoires=analyzed_repertoires,
                                                 statDistance=WassersteinStatDistance(combo[0], combo[1])
                                                )
            study.optimize(func=objectiveFunction,n_trials=3)
            results[repertoire_group] = study

            # Best result
            print("Best score:", study.best_value)
            print("Best params:", study.best_params)
            print("Generating sample networks") 
            for repertoire_dataset in all_repertoires[repertoire_group]:
                sampleRepertoire = all_repertoires[repertoire_group][repertoire_dataset][0]
                immuneNet = makeBestNetwork(sampleRepertoire, study)
                ImmuneNetworkMapper.toPickle(network=immuneNet,path="network_{combo[0]}_{combo[1]}_{repertoire_group}_{repertoire_dataset}_{runId}.csv")

        with open(f"results_{combo[0]}_{combo[1]}_{runId}", "wb") as f:
            pickle.dump(results, f)






