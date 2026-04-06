import optuna
from pathlib import Path
from itertools import combinations
import numpy as np
import uuid
import pickle
import argparse
import os
import random
from multiprocessing import Process, Pool;

from immune_nets.analysis.optunaObjectives.initialAnalysisOjectiveBuilder import initialAnalysisObjectiveBuilder 

# from immune_nets.models import Repertoire,Dataset
# from immune_nets.db import engine
# from sqlalchemy.orm import Session, selectinload
# from sqlalchemy import insert,select,delete
from immune_nets.analysis.statDistances import WassersteinStatDistance, PairwiseDistributionDistance
from immune_nets.analysis.scoringParadigms import PairwiseScoringParadigm
from immune_nets.mappers import ImmuneNetworkMapper
from immune_nets.factories import ImmuneRepertoireFactory

from immune_nets.entities import ImmuneRepertoire
from immune_nets.entities import ImmuneNetwork
from immune_nets.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from immune_nets.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from immune_nets.creation.distance.levenshtein import levenshteinDistance
from immune_nets.creation.distance.alignment import sequenceAligner
from immune_nets.cluster.utils.commonMethods import loadDatasetsFromFile, createTestGroups

groups = ["covid","healthy"]
groupSize = [2,2]



def stopIfThresholdReached(study, trial):
    if trial.value is not None and trial.value >= 1.0:
        print("Max score reached - stopping heuristic")
        study.stop()


        
def runClusterJob(repertoireTestGroup, test_group, distributionName, run_id, rank, numOfTrials=100, distributionComputingPoolSize=3, testCase=False):
    distributionNames = ["degreeDistribution", "componentSizeDistribution", "componentProportionDistribution"]

    print(f"Process {rank} is starting computation for {distributionName}.")
    results = {}
    
    # print(f"Process {rank} is running study for {test_group} repertoires")

    analyzed_repertoires = repertoireTestGroup
    distanceType = PairwiseDistributionDistance(distributionName)
    scoringParadigm = PairwiseScoringParadigm(distanceType)
    objectiveFunction = initialAnalysisObjectiveBuilder(repertoires=analyzed_repertoires,
                                         scoringParadim=scoringParadigm,
                                         rank=rank
                                        )


    study = optuna.create_study(direction="maximize")
    study.optimize(func=objectiveFunction, n_trials=numOfTrials, callbacks=[stopIfThresholdReached])
    print("Process {world_rank}: Best score:", study.best_value)
    print("Process {world_rank}: Best params:", study.best_params)

    # Best result

    filename = f"results_initial_analysis_{distributionName}_{test_group}_{run_id}.pkl"
    if not testCase: 
        with open(filename, "wb") as f:
            pickle.dump(results, f)
    else:
        rand = random.randint(0, 1_000_000)
        print(f"This is testcase. Process rank is {rank}. Id is {rand} and thus filename is {filename}.")

    print(f"Finished processing for {distributionName}")



def runProcesses(repertoireDatasets, run_id, numOfTrials=100, testCase=False):
    distributionNames = ["degreeDistribution", "componentSizeDistribution", "componentProportionDistribution"]
    processes = []
    repertoireTestGroups = createTestGroups(repertoireDatasets)
   
    for rank, distributionName in enumerate(distributionNames):
        for test_group in repertoireTestGroups:
            p = Process(target=runClusterJob,
                        args=(repertoireTestGroups[test_group], test_group, distributionName, run_id, rank, numOfTrials, testCase))
            p.start()
            processes.append(p)

    for p in processes:
        p.join()




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("--group-paths", type=str, required=True)
    parser.add_argument("--run-id", type=int, required=True)
    args = parser.parse_args()

    run_id = args.run_id

    groupPaths = args.group_paths
    groupPaths = groupPaths.split(",")
    groupPaths = { pair.split(":")[0]:pair.split(":")[1] for pair in groupPaths}

    all_repertoires = loadDatasetsFromFile(groupPaths)

    runProcesses(all_repertoires, run_id)









