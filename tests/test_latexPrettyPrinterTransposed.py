import pickle
import re
import os

from src.presentation import LatexPrettyPrinterTransposed

from src.analysis.visualization.multiGraphChart import multiGraphChart
from tqdm import tqdm

from pathlib import Path

def find_project_root(start: Path = Path(__file__)):
    for p in [start, *start.parents]:
        if (p / ".git").exists():
            return p
    raise RuntimeError("Project root not found")

home = os.environ.get("HOME")

dataDir = f"{home}/Documents/studia/heuristicResults/" 

rootPath = find_project_root()


resultList = [
    "results_componentProportionDistribution_6e3e3a9d-4abf-445a-b887-d5eaf9ee15af",
    "results_componentSizeDistribution_e7031d3d-c12d-454e-9c36-2d856eb4a843",
    "results_degreeDistribution_3e402735-6b9c-492f-abdb-6205f1caca28",
]



optResultFilenames = [
    "results_optimizer_componentProportionDistribution_1.pkl",
    "results_optimizer_componentProportionDistribution_3.pkl",
    "results_optimizer_componentProportionDistribution_4.pkl",
    "results_optimizer_componentProportionDistribution_5.pkl",
    "results_optimizer_componentProportionDistribution_6.pkl",
    "results_optimizer_componentSizeDistribution_1.pkl",
    "results_optimizer_componentSizeDistribution_3.pkl",
    "results_optimizer_componentSizeDistribution_4.pkl",
    "results_optimizer_componentSizeDistribution_5.pkl",
    "results_optimizer_componentSizeDistribution_6.pkl",
    "results_optimizer_degreeDistribution_1.pkl",
    "results_optimizer_degreeDistribution_3.pkl",
    "results_optimizer_degreeDistribution_4.pkl",
    "results_optimizer_degreeDistribution_5.pkl",
    "results_optimizer_degreeDistribution_6.pkl",
]

distributionNameDict = {
    "degreeDistribution": "rozkładu stopni wierzchołków",
    "componentSizeDistribution": "rozkładu rozmiaru komponentów",
    "componentProportionDistribution": "rozkładu proporcji komponentów",
}

class BestTrial:
    def __init__(self, number,  value, params, user_attrs):
        self.number = number
        self.values = value
        self.params = params 
        self.user_attrs = user_attrs

optResultDict = { distributionName:[] for distributionName in distributionNameDict}

for resultFilename in optResultFilenames:
    resultFilenameParts = resultFilename.split("_")
    distributionName = resultFilenameParts[2]
    with open(rootPath / resultFilename, "rb") as f:
        result = pickle.load(f)
        isCustom = False if not int(resultFilenameParts[3][0]) in [1,3] else True
        bestTrial = result if isCustom else result.best_trial
        optResultDict[distributionName].append(bestTrial)
        
    


networkPaths = [
    'network_componentProportionDistribution_covid_1 vs covid_2_covid_1_6e3e3a9d-4abf-445a-b887-d5eaf9ee15af.csv',
    'network_componentProportionDistribution_covid_1 vs covid_2_covid_2_6e3e3a9d-4abf-445a-b887-d5eaf9ee15af.csv',
    'network_componentProportionDistribution_healthy_1 vs healthy_2_healthy_1_6e3e3a9d-4abf-445a-b887-d5eaf9ee15af.csv',
    'network_componentProportionDistribution_healthy_1 vs healthy_2_healthy_2_6e3e3a9d-4abf-445a-b887-d5eaf9ee15af.csv',
    'network_componentProportionDistribution_healthy vs covid_covid_6e3e3a9d-4abf-445a-b887-d5eaf9ee15af.csv',
    'network_componentProportionDistribution_healthy vs covid_healthy_6e3e3a9d-4abf-445a-b887-d5eaf9ee15af.csv',
    'network_componentSizeDistribution_covid_1 vs covid_2_covid_1_e7031d3d-c12d-454e-9c36-2d856eb4a843.csv',
    'network_componentSizeDistribution_covid_1 vs covid_2_covid_2_e7031d3d-c12d-454e-9c36-2d856eb4a843.csv',
    'network_componentSizeDistribution_healthy_1 vs healthy_2_healthy_1_e7031d3d-c12d-454e-9c36-2d856eb4a843.csv',
    'network_componentSizeDistribution_healthy_1 vs healthy_2_healthy_2_e7031d3d-c12d-454e-9c36-2d856eb4a843.csv',
    'network_componentSizeDistribution_healthy vs covid_covid_e7031d3d-c12d-454e-9c36-2d856eb4a843.csv',
    'network_componentSizeDistribution_healthy vs covid_healthy_e7031d3d-c12d-454e-9c36-2d856eb4a843.csv',
    'network_degreeDistribution_covid_1 vs covid_2_covid_1_3e402735-6b9c-492f-abdb-6205f1caca28.csv',
    'network_degreeDistribution_covid_1 vs covid_2_covid_2_3e402735-6b9c-492f-abdb-6205f1caca28.csv',
    'network_degreeDistribution_healthy_1 vs healthy_2_healthy_1_3e402735-6b9c-492f-abdb-6205f1caca28.csv',
    'network_degreeDistribution_healthy_1 vs healthy_2_healthy_2_3e402735-6b9c-492f-abdb-6205f1caca28.csv',
    'network_degreeDistribution_healthy vs covid_covid_3e402735-6b9c-492f-abdb-6205f1caca28.csv',
    'network_degreeDistribution_healthy vs covid_healthy_3e402735-6b9c-492f-abdb-6205f1caca28.csv',
]

allNetworks = {}

groups = [
    'healthy vs covid',
    'covid_1 vs covid_2',
    'healthy_1 vs healthy_2'
]

for networkFilename in networkPaths:
    networkPath = dataDir + networkFilename
    splitFilename = networkFilename.split("_")
    distributionName = splitFilename[1]
    group = "dummy"
    for checked_group in groups:
        if re.search(checked_group, networkFilename):
            group = checked_group
            break

    networkGroup = splitFilename[-2] if not splitFilename[-2] in ["1","2"] else (splitFilename[-3] + "_" + splitFilename[-2])

    if not distributionName in allNetworks:
        allNetworks[distributionName] = {}
    if not group in allNetworks[distributionName]:
        allNetworks[distributionName][group] = {} 
    with open(networkPath, "rb") as f:
        allNetworks[distributionName][group][networkGroup] = pickle.load(f)

# print(allNetworks)



def generateBasicTables():

    print(f"\\section{{Analiza rozkładów związanych z generowaniem sieci}}")
    for resultFilename in resultList:
        resultPath = dataDir + resultFilename
        with open(resultPath, "rb") as f:
            result = pickle.load(f)
        result = { resname:result[resname] for resname in result if resname != "universal"}
        distributionName = resultFilename.split("_")[1]

        LatexPrettyPrinterTransposed().printTable(result, allNetworks[distributionName], distributionNameDict[distributionName], distributionName)
        # print(LatexPrettyPrinterTransposed().generateLatexGraphHeader(result))

def generateBasicNetworks():
    for distributionName in allNetworks:
        print(f"Calculating graphs for {distributionName}")
        for test_group in tqdm(allNetworks[distributionName]):
            networksDict = allNetworks[distributionName][test_group]
            # plotTitles = [ f"Sieć próbki {group} dla {distributionNameDict[distributionName]}" for group in networksDict] 
            plotTitles = ["A", "B"]
            networks = [ networksDict[group] for group in networksDict ]
            multiGraphChart(plotTitles, networks, f"{distributionName}_{test_group}.png") 



def saveResults():
    num = 3
    distribution = "componentProportionDistribution"
    number = 56
    value = [0.5578741867428201]
    params = {'threshold': 0.392900516710026, 'distance': 'alignment', 'algorithm_name': 'simpleBetaDistance', 'substitution_matrix': 'PAM250'}
    user_attrs = {'healthy vs covid score': 0.339665861432621 , 'covid_1 vs covid_2 score': 0.11284417096583933, 'healthy_1 vs healthy_2 score': 0.18750591947735692, 'between': 0.339665861432621, 'within': 0.150175045221598}
    bestTrial = BestTrial(number, value, params, user_attrs)
    with open(f"results_optimizer_{distribution}_{num}.pkl", "wb") as f:
        pickle.dump(bestTrial, f)

    distribution = "componentSizeDistribution"
    number = 0
    value = [0.0]
    params = {'threshold': 0.2500209182377269, 'distance': 'levenshtein', 'algorithm_name': 'simpleBetaDistance'}
    user_attrs = {'healthy vs covid score': 0.0 , 'covid_1 vs covid_2 score': 0.0, 'healthy_1 vs healthy_2 score': 0.0, 'between': 0.0, 'within': 0.6}
    bestTrial = BestTrial(number, value, params, user_attrs)
    with open(f"results_optimizer_{distribution}_{num}.pkl", "wb") as f:
        pickle.dump(bestTrial, f)


    distribution = "degreeDistribution"
    number = 33
    value = [0.5145577452426411]
    params = {'threshold': 0.2014423281575187, 'distance': 'alignment', 'algorithm_name': 'simpleVectorBetaDistance', 'substitution_matrix': 'BLOSUM62'}
    user_attrs = {'healthy vs covid score': 0.6 , 'covid_1 vs covid_2 score': 0.5825307057088307, 'healthy_1 vs healthy_2 score': 0.0, 'between': 0.6, 'within': 0.291265352854415}
    bestTrial = BestTrial(number, value, params, user_attrs)
    with open(f"results_optimizer_{distribution}_{num}.pkl", "wb") as f:
        pickle.dump(bestTrial, f)

    num = 1
    distribution = "componentProportionDistribution"
    number = 75
    value = [0.618408267991226]
    params = {'threshold': 0.38573678908300957, 'distance': 'alignment', 'algorithm_name': 'simpleBetaDistance', 'substitution_matrix': 'PAM250'}
    user_attrs = {'healthy vs covid score': 0.334458822155024, 'covid_1 vs covid_2 score': 0.0964205635852593, 'healthy_1 vs healthy_2 score': 0.158832878878241, 'between': 0.334458822155024, 'within': 0.12762672123175}
    bestTrial = BestTrial(number, value, params, user_attrs)
    with open(f"results_optimizer_{distribution}_{num}.pkl", "wb") as f:
        pickle.dump(bestTrial, f)

    distribution = "componentSizeDistribution"
    number = 0
    value = [0.0]
    params = {'threshold': 0.27577915272646125, 'distance': 'alignment', 'algorithm_name': 'simpleBetaDistance', 'substitution_matrix': 'PAM30'}
    user_attrs = {'healthy vs covid score': 0.0 , 'covid_1 vs covid_2 score': 0.0, 'healthy_1 vs healthy_2 score': 0.0, 'between': 0.0, 'within': 0.6}
    bestTrial = BestTrial(number, value, params, user_attrs)
    with open(f"results_optimizer_{distribution}_{num}.pkl", "wb") as f:
        pickle.dump(bestTrial, f)

    distribution = "degreeDistribution"
    number = 94
    value = [0.653506187165161]
    params = {'threshold': 0.2014423281575187, 'distance': 'alignment', 'algorithm_name': 'simpleVectorBetaDistance', 'substitution_matrix': 'BLOSUM62'}
    user_attrs = {'healthy vs covid score': 0.7 , 'covid_1 vs covid_2 score': 0.4850913379687746, 'healthy_1 vs healthy_2 score': 0.0, 'between': 0.7, 'within': 0.242545668984387}
    bestTrial = BestTrial(number, value, params, user_attrs)
    with open(f"results_optimizer_{distribution}_{num}.pkl", "wb") as f:
        pickle.dump(bestTrial, f)

def generateResultTables():
    print(f"\\section{{Porównanie parametrów generowanych na podstawie różnych rozkładów}}")
    for distributionName in optResultDict:
        result = optResultDict[distributionName]
        LatexPrettyPrinterTransposed().printOptTable(result, distributionNameDict[distributionName])





import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables", action="store_true")
    parser.add_argument("--networks", action="store_true")
    parser.add_argument("--save-results", action="store_true")
    parser.add_argument("--result-tables", action="store_true")
    args = parser.parse_args()
    if args.tables:
        generateBasicTables()
    if args.networks:
        generateBasicNetworks()
    if args.save_results:
        saveResults()
    if args.result_tables:
        generateResultTables()
        


