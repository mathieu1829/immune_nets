import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA
import umap
from sklearn.preprocessing import MinMaxScaler
import numpy as np

from src.analysis.visualization.graphVisualization import graphVisualization
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.immuneRepertoire import ImmuneRepertoire
from src.creation.distance.alignment import sequenceAligner
from src.analysis.methods.graphletComposition import graphletComposition
from sklearn.preprocessing import MinMaxScaler

def graphStatisticBarplots(immuneNets):
    statList = [
                "vertice_num",
                "isolated_vertices_num",
                "edge_density",
                "percolation_threshold",
                "density",
                "eccentrity.mean()",
                "eigenvector_centrality.mean()",
                "harmonic_centrality.mean()",
                "giant_component",
                "betweenness.mean()",
                "diameter",
                "mean_closeness",
                "mean_shortest_path",
                "expected_pagerank",
                "expected_degree",
                "component_count",
                "expected_component_size"
            ]
    immuneNetsStats = { group:graphletComposition(immuneNets[group]).toList() for group in immuneNets}
    
    groups = [ group for group in immuneNets ]
    x = np.arange(3)
    width = 0.25
    
    for stat_num, stat in enumerate(statList) :
        currentStat = [ immuneNetsStats[group][stat_num] for group in immuneNetsStats]
        plt.bar(x+stat_num*width, currentStat, width, label=stat)
    plt.xticks(x + width, labels = groups) 

    plt.show()

           
if __name__ == "__main__":
    groups = ["leukemia", "covid", "healthy"]

    root_dir = Path(__file__).parent.parent.parent.parent

    leukemia_path = root_dir  / "tests/test_data/leukemia_test_clonotypes.csv" # leukemia
    covid_path = root_dir / "tests/test_data/covid_test_clonotypes.csv" # covid
    healthy_path = root_dir / "tests/test_data/healthy_test_clonotypes_1.csv" #healthy

    # TO DO - get repertoires for each group from database
    repertoire_list = [
            ImmuneRepertoire.fromCSV(path=leukemia_path, name="leukemia", desc=" "),
            ImmuneRepertoire.fromCSV(path=covid_path, name="leukemia", desc=" "),
            ImmuneRepertoire.fromCSV(path=healthy_path, name="leukemia", desc=" "),
            ]
    distance_fun = sequenceAligner("BLOSUM62")
    grouped_immuneNets = { 
                   group: simpleBetaDistance(repertoire=repertoire_list[i],
                                                distance=distance_fun,
                                                threshold=0.2)
                         
                   for i,group in enumerate(groups)
                   }
    graphStatisticBarplots(grouped_immuneNets)

