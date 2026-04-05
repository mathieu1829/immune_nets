import matplotlib.pyplot as plt
from pathlib import Path

from immune_nets.analysis.visualization.graphVisualization import graphVisualization
from immune_nets.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from immune_nets.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from immune_nets.creation.distance.alignment import sequenceAligner
from immune_nets.factories import ImmuneRepertoireFactory


def singleGraphChart(name, immuneNet, filename):
    
    fig, ax = plt.subplots(figsize=(16, 16))

    ax.set_title(name)
    graphVisualization( immuneNet, ax)

    plt.savefig(filename)

if __name__ == "__main__":
    root_dir = Path(__file__).parent.parent.parent.parent
    
    healthy_path = root_dir / "tests/test_data/healthy_test_clonotypes_1.csv" #healthy
    repertoire = ImmuneRepertoireFactory.fromCSVTest(healthy_path)
    distance_fun = sequenceAligner("BLOSUM62")
    immuneNet = simpleVectorBetaDistance(repertoire=repertoire,
                                     distance=distance_fun,
                                     threshold=0.2)

    singleGraphChart("healthy repertorie",repertoire,immuneNet)

