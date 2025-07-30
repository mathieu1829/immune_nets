import matplotlib.pyplot as plt
from pathlib import Path

from src.analysis.visualization.graphVisualization import graphVisualization
from src.creation.io_strategies.test_csv_strategy import *
from src.creation.algorithms.simple_beta_distance import simple_beta_distance
from src.creation.distance.alignment import sequenceAligner


def singleGraphChart(name, repertoire, immuneNet):
    
    fig, ax = plt.subplots(figsize=(16, 16))

    ax.set_title(name)
    graphVisualization(repertoire, immuneNet, ax)

    plt.show()

if __name__ == "__main__":
    root_dir = Path(__file__).parent.parent.parent.parent
    
    healthy_path = root_dir / "tests/test_data/healthy_test_clonotypes_1.csv" #healthy
    repertoire = test_csv_strategy().input(healthy_path)
    distance_fun = sequenceAligner("BLOSUM62")
    immuneNet = simple_beta_distance(repertoire=repertoire,
                                     distance=distance_fun,
                                     threshold=0.2)

    singleGraphChart("healthy repertorie",repertoire,immuneNet)

