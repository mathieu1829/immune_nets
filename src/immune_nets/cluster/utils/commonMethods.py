import os
from immune_nets.factories import ImmuneRepertoireFactory
from immune_nets.factories import RepertoireFactory

from immune_nets.entities import ImmuneRepertoire
from immune_nets.entities import ImmuneNetwork
from immune_nets.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from immune_nets.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from immune_nets.creation.distance.levenshtein import levenshteinDistance
from immune_nets.creation.distance.alignment import sequenceAligner



def divideIntoSubgroups(groupedRepertoires, group):
  repertoireGroup = groupedRepertoires[group]
  mid = len(repertoireGroup) // 2
  return {f"{group}_1":repertoireGroup[:mid], f"{group}_2":repertoireGroup[mid:]}

def loadDatasetsFromFile(groupPaths, db=False):

    repertoireDatasets = {}

    for group in groupPaths:
        repertoireList = []
        for file in os.listdir(groupPaths[group]):
            path = groupPaths[group] + "/" + file
            metadaGroups = ["group", "id", "description"]
            metadata = { group:data for group, data in zip(metadaGroups, file.split("_"))}
            args = {"name":f"{metadata['group']} {metadata['id']}",
                    "desc":f"{metadata['description']}",
                    "path":path
                    }
            
            repertoire = ImmuneRepertoireFactory.fromCSV(**args) if not db else RepertoireFactory.fromCSV(**args)
            repertoireList.append(repertoire)
        repertoireDatasets[group] = repertoireList
    return repertoireDatasets

def createTestGroups(repertoireStatsDatasets):
    groups = [group for group in repertoireStatsDatasets]
    testGroups = { f"healthy vs {group}":{"healthy":repertoireStatsDatasets["healthy"], group:repertoireStatsDatasets[group]} for group in groups if not group == "healthy"}
    # all_repertoires["universal"] = dataset_repertoires
    for group in groups:
      testGroups[f"{group}_1 vs {group}_2"] = divideIntoSubgroups(repertoireStatsDatasets, group)
    return testGroups

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
