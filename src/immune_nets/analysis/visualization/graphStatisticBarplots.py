import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA
import umap
from sklearn.preprocessing import MinMaxScaler
import numpy as np

from immune_nets.analysis.visualization.graphVisualization import graphVisualization
from immune_nets.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from immune_nets.entities import ImmuneRepertoire
from immune_nets.creation.distance.alignment import sequenceAligner
from immune_nets.entities import GraphStats
from immune_nets.mappers import GraphStatsMapper
from sklearn.preprocessing import MinMaxScaler
from sqlalchemy.orm import Session, selectinload
from sqlalchemy import insert,select,delete
from immune_nets.db import engine
from immune_nets.models import Network
from immune_nets.mappers import NetworkMapper

def graphStatisticBarplots(immuneNets):
    immuneNetsStats = { group:GraphStatsMapper.toStatVector(GraphStats(immuneNets[group])) for group in immuneNets}

    stds = []
    means = []
    for i,stat in enumerate(GraphStatsMapper.colnames()):
        statCol = [immuneNetsStats[group][i] for group in immuneNetsStats]
        scaler = MinMaxScaler()
        scaledCol = scaler.fit_transform([[v] for v in statCol])
        scaledCol = [float(x[0]) for x in scaledCol]
        stds.append(np.std(scaledCol))
        means.append(np.mean(scaledCol))

    x = np.arange(len(stds))
    # width = 0.05                      

    for name,stat in zip(["std","mean"],[stds,means]):
        plt.bar(x, stat, color='skyblue')

        plt.xticks(x, GraphStatsMapper.colnames(), rotation=45, ha='right')

        plt.ylabel("Value")
        plt.title(f"{name} Barplot")

        for i, v in enumerate(stat):
            plt.text(x[i], v + 0.01, f"{v:.2f}", ha='center', va='bottom')

        plt.tight_layout()
        plt.show()

    # for i,stat in enumerate(statList):
    #     statCol = [immuneNetsStats[group][i] for group in immuneNetsStats]
    #     scaler = MinMaxScaler()
    #     scaledCol = scaler.fit_transform([[v] for v in statCol])
    #     scaledCol = [float(x[0]) for x in scaledCol]
    #     for idx,group in enumerate(immuneNetsStats):
    #         immuneNetsStats[group][i] = scaledCol[idx]
    # 
    # x = np.arange(len(groups))        
    # width = 0.05                      
    #
    # plt.figure(figsize=(14, 6))
    #
    # for i, stat in enumerate(statList):
    #     stat_values = [immuneNetsStats[group][i] for group in groups]
    #     plt.bar(x + i * width - (len(statList) / 2) * width, stat_values, width, label=stat)
    #
    # plt.xticks(x, groups)
    # plt.ylabel("Scaled statistic (MinMax)")
    # plt.title("Network Statistics per Group")
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2)
    # plt.tight_layout()
    # plt.show()

           
if __name__ == "__main__":
    groups = ["leukemia", "covid", "healthy"]

    grouped_immuneNets = {}
    with Session(engine) as session: 
        for group in groups :
            stmt = select(Network).options(selectinload(Network.network_edges)).where(Network.name == f"{group} test network")
            result = session.execute(stmt)
            network = result.scalars().first()
            grouped_immuneNets[group] = NetworkMapper.toImmuneNetwork(network)
            print(grouped_immuneNets[group].graph)


    graphStatisticBarplots(grouped_immuneNets)

