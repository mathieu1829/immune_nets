import optuna
from pathlib import Path
from itertools import combinations
import numpy as np
import uuid
import pickle
import argparse
import os
from multiprocessing import Process, Pool;
import random

from src.analysis.optunaObjectives import optimizerObjectiveBuilder

# from src.models import Repertoire,Dataset
# from src.db import engine
# from sqlalchemy.orm import Session, selectinload
# from sqlalchemy import insert,select,delete
from src.analysis.statDistances import WassersteinStatDistance, PairwiseDistributionDistance
from src.analysis.scoringParadigms import PairwiseScoringParadigm
from src.mappers import ImmuneNetworkMapper
from src.factories import ImmuneRepertoireFactory

from src.entities import ImmuneRepertoire
from src.entities import ImmuneNetwork
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from src.creation.distance.levenshtein import levenshteinDistance
from src.creation.distance.alignment import sequenceAligner
from src.analysis.optunaObjectives import compareGroups


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

# dataset_repertoires = {}
# with Session(engine) as session: 
#     for group in groups :
#         stmt = select(Dataset).options(selectinload(Dataset.repertoires).selectinload(Repertoire.clonotypes)).where(Dataset.name == f"{group} test dataset")
#         result = session.execute(stmt)
#         dataset = result.scalars().first()
#         dataset_repertoires[group] = [ RepertoireMapper.toImmuneRepertoire(repertoire) for repertoire in dataset.repertoires]
#
# all_repertoires = { f"healthy vs {group}":{"healthy":dataset_repertoires["healthy"], group:dataset_repertoires[group]} for group in groups if not group == "healthy"}
# all_repertoires["universal"] = dataset_repertoires
# for group in groups:
#   all_repertoires[f"{group}_1 vs {group}_2"] = shuffleGroupAndDivide(dataset_repertoires, group)

def loadClusterJobDatasetsFromFile(groupPaths):

    dataset_repertoires = {}

    for group in groupPaths:
        repertoireList = []
        for file in os.listdir(groupPaths[group]):
            path = groupPaths[group] + "/" + file
            metadaGroups = ["group", "id", "description"]
            metadata = { group:data for group, data in zip(metadaGroups, file.split("_"))}
            repertoireList.append(ImmuneRepertoireFactory.fromCSV(name=f"{metadata['group']} {metadata['id']}",desc=f"{metadata['description']}",path=path))
        dataset_repertoires[group] = repertoireList

    all_repertoires = { f"healthy vs {group}":{"healthy":dataset_repertoires["healthy"], group:dataset_repertoires[group]} for group in groups if not group == "healthy"}
    # all_repertoires["universal"] = dataset_repertoires
    for group in groups:
      all_repertoires[f"{group}_1 vs {group}_2"] = shuffleGroupAndDivide(dataset_repertoires, group)
    return all_repertoires

def stopIfThresholdReached(study, trial):
    if trial.value is not None and trial.value >= 1.0:
        print("Max score reached - stopping heuristic")
        study.stop()

        
def runClusterJob(allRepertoires, distributionName, run_id, rank, numOfTrials=100, testCase=False):
    
    print(f"Process {rank} is starting computation for {distributionName}.")
    
    distanceType = PairwiseDistributionDistance(distributionName)
    scoringParadigm = PairwiseScoringParadigm(distanceType)

    study = optuna.create_study(direction="maximize")
    objectiveFunction = optimizerObjectiveBuilder(allRepertoires=allRepertoires,
                                         scoringParadigm=scoringParadigm,
                                         rank=rank
                                        )
    study.optimize(func=objectiveFunction,n_trials=numOfTrials, callbacks=[stopIfThresholdReached])

    # Best result
    print(f"Process {rank}: Best score: {study.best_value}")
    print(f"Process {rank}: Best params: {study.best_params}" )

    if not testCase:  
        with open(f"results_optimizer_{distributionName}_{run_id}.pkl", "wb") as f:
            pickle.dump(study, f)
    else:
        filename = f"proc_{rank}_{run_id}.txt"
        print(f"This is testcase. Process rank is {rank}. Id is {run_id} and thus filename is {filename}.")

        with open(filename, "w") as f:
            f.write("a")

    print(f"Finished processing for {distributionName}")



        
    



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("--group-paths", type=str, required=True)
    parser.add_argument("--run-id", type=int, required=True)
    args = parser.parse_args()

    run_id = args.run_id

    groupPaths = args.group_paths
    groupPaths = groupPaths.split(",")
    groupPaths = { pair.split(":")[0]:pair.split(":")[1] for pair in groupPaths}

    all_repertoires = loadClusterJobDatasetsFromFile(groupPaths)

    distributionNames = ["degreeDistribution", "componentSizeDistribution", "componentProportionDistribution"]
    processes = []
   
    for distributionName in distributionNames:
        p = Process(target=runClusterJob,
                    args=(all_repertoires, distributionName, run_id))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()








