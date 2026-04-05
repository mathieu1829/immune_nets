import pickle
import sys
import re
import os

sys.path.append("../")

from immune_nets.presentation import LatexPrettyPrinter

home = os.environ.get("HOME")
print(home)

dataDir = f"{home}/Documents/studia/heuristicResults/" 

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
    'network_proportionCountDistribution_covid_1 vs covid_2_covid_1_53a87b52-22d8-4850-956d-67beb45c7f17.csv',
    'network_proportionCountDistribution_covid_1 vs covid_2_covid_2_53a87b52-22d8-4850-956d-67beb45c7f17.csv',
    'network_proportionCountDistribution_healthy_1 vs healthy_2_healthy_1_53a87b52-22d8-4850-956d-67beb45c7f17.csv',
    'network_proportionCountDistribution_healthy_1 vs healthy_2_healthy_2_53a87b52-22d8-4850-956d-67beb45c7f17.csv',
    'network_proportionCountDistribution_healthy vs covid_covid_53a87b52-22d8-4850-956d-67beb45c7f17.csv',
    'network_proportionCountDistribution_healthy vs covid_healthy_53a87b52-22d8-4850-956d-67beb45c7f17.csv',
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

distributionNameDict = {
    "degreeDistribution": "rozkładu stopni wierzchołków",
    "componentSizeDistribution": "rozkładu rozmiaru komponentów",
    "proportionCountDistribution": "rozkładu proporcji wierzchołków",
    "componentProportionDistribution": "rozkładu proporcji komponentów",
 
}

from immune_nets.analysis.visualization.multiGraphChart import multiGraphChart
from tqdm import tqdm

for distributionName in allNetworks:
    print(f"Calculating graphs for {distributionName}")
    for test_group in tqdm(allNetworks[distributionName]):
        networksDict = allNetworks[distributionName][test_group]
        plotTitles = [ f"Sieć próbki {group} dla {distributionNameDict[distributionName]}" for group in networksDict] 
        networks = [ networksDict[group] for group in networksDict ]
        multiGraphChart(plotTitles, networks, f"{distributionName}_{test_group}.png") 
