import optuna
from pathlib import Path
from itertools import combinations
import numpy as np
import pickle
import argparse
from multiprocessing import Process, Pool;

from immune_nets.analysis.optunaObjectives import optimizerObjectiveBuilder

from immune_nets.analysis.statDistances import WassersteinStatDistance, PairwiseDistributionDistance
from immune_nets.analysis.scoringParadigms import PairwiseScoringParadigm
from immune_nets.cluster.utils.commonMethods import loadDatasetsFromFile
from immune_nets.mappers import ImmuneNetworkMapper

groups = ["covid","healthy"]
groupSize = [2,2]



def stopIfThresholdReached(study, trial):
    if trial.value is not None and trial.value >= 1.0:
        print("Max score reached - stopping heuristic")
        study.stop()

        
def runClusterJob(allRepertoires, distributionName, run_id, rank, numOfTrials=100, testCase=False):
    
    print(f"Process {rank} is starting computation for {distributionName}.")
    
    distanceType = PairwiseDistributionDistance(distributionName)
    scoringParadigm = PairwiseScoringParadigm(distanceType)

    study = optuna.create_study(direction="maximize")
    objectiveFunction = optimizerObjectiveBuilder(repertoireDatasets=allRepertoires,
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



def runProcesses(repertoireDatasets, run_id, numOfTrials=100, testCase=False):
    distributionNames = ["degreeDistribution", "componentSizeDistribution", "componentProportionDistribution"]
    processes = []
   
    for rank, distributionName in enumerate(distributionNames):
        p = Process(target=runClusterJob,
                    args=(repertoireDatasets, distributionName, run_id, rank, numOfTrials, testCase))
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

    repertoireDatasets = loadDatasetsFromFile(groupPaths)

    runProcesses(repertoireDatasets, run_id)










