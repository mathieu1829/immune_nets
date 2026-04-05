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
from immune_nets.factories import ImmuneRepertoireFactory

def graphSimilarityChart(grouped_immuneNets):
    group_results = []
    group_sizes = [len(grouped_immuneNets[group]) for group in grouped_immuneNets]
    total_samples = sum(group_sizes)
    for i in range(1,len(group_sizes)):
        group_sizes[i] = group_sizes[i] + group_sizes[i-1]

    for group in grouped_immuneNets: 
        group_immuneNets = grouped_immuneNets[group]

        for immuneNet in group_immuneNets:
            stats = GraphStats(immuneNet)
            group_results.append(GraphStatsMapper.toList(stats))
    scaler = MinMaxScaler()
    normalizedResults = scaler.fit_transform(group_results)
    if total_samples > 3:
        n_neighbors = max(2, min(15, total_samples - 1))  # safe, auto-adjust
        reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=0.1, metric='cosine')
        reduced = reducer.fit_transform(np.array(normalizedResults))
        title = "UMAP"
    else:
        reduced = PCA(n_components=2).fit_transform(np.array(normalizedResults))
        title = "PCA"



    for i,group in enumerate(grouped_immuneNets):
        start = group_sizes[i-1] if i > 0 else 0
        end = group_sizes[i]
        plt.scatter(reduced[start:end, 0], reduced[start:end, 1], alpha=0.6, label=group)
    
    plt.title(f"2D Graph Embedding Visualization ({title})")
    plt.legend(title="Group")
    plt.xlabel(f"{title}-1")
    plt.ylabel(f"{title}-2")
    plt.show()
            
if __name__ == "__main__":
    groups = ["leukemia", "covid", "healthy"]

    root_dir = Path(__file__).parent.parent.parent.parent

    leukemia_path = root_dir  / "tests/test_data/leukemia_test_clonotypes_0.csv" # leukemia
    covid_path = root_dir / "tests/test_data/covid_test_clonotypes_0.csv" # covid
    healthy_path = root_dir / "tests/test_data/healthy_test_clonotypes_1.csv" #healthy

    # TO DO - get repertoires for each group from database
    repertoire_list = [
            ImmuneRepertoireFactory.fromCSVTest(leukemia_path),
            ImmuneRepertoireFactory.fromCSVTest(covid_path),
            ImmuneRepertoireFactory.fromCSVTest(healthy_path),
            ]
    distance_fun = sequenceAligner("BLOSUM62")
    grouped_immuneNets = { 
                   group:[
                           simpleBetaDistance(repertoire=repertoire_list[i],
                                                distance=distance_fun,
                                                threshold=0.2)
                       ]  
                   for i,group in enumerate(groups)
                   }
    graphSimilarityChart(grouped_immuneNets)

