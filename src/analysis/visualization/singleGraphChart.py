import matplotlib.pyplot as plt
from pathlib import Path

from src.analysis.visualization.graphVisualization import graphVisualization
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from src.creation.distance.alignment import sequenceAligner
from src.factories import ImmuneRepertoireFactory


def singleGraphChart(name, repertoire, immuneNet):
    
    fig, ax = plt.subplots(figsize=(16, 16))

    ax.set_title(name)
    graphVisualization(repertoire, immuneNet, ax)

    plt.show()

if __name__ == "__main__":
    root_dir = Path(__file__).parent.parent.parent.parent
    
    healthy_path = root_dir / "tests/test_data/healthy_test_clonotypes_1.csv" #healthy
    repertoire = ImmuneRepertoireFactory.fromCSVTest(healthy_path)
    distance_fun = sequenceAligner("BLOSUM62")
    immuneNet = simpleVectorBetaDistance(repertoire=repertoire,
                                     distance=distance_fun,
                                     threshold=0.2)

    singleGraphChart("healthy repertorie",repertoire,immuneNet)

