import matplotlib.pyplot as plt
from pathlib import Path

from src.analysis.visualization.graphVisualization import graphVisualization
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.distance.alignment import sequenceAligner
from src.factories import ImmuneRepertoireFactory

def multiGraphChart(groups, repertoires, immuneNets):
    fig, axes = plt.subplots(nrows=1,ncols=len(groups),figsize=(16, 16))

    for group_num, group in enumerate(groups):
        axes[group_num].set_title(group)
        graphVisualization(repertoires[group_num],
                           immuneNets[group_num],
                           axes[group_num]
                           )

    plt.show()


if __name__ == "__main__":

    groups = ["leukemia", "covid", "healthy"]
    root_dir = Path(__file__).parent.parent.parent.parent
    
    leukemia_path = root_dir  / "tests/test_data/leukemia_test_clonotypes.csv" # leukemia
    covid_path = root_dir / "tests/test_data/covid_test_clonotypes.csv" # covid
    healthy_path = root_dir / "tests/test_data/healthy_test_clonotypes_1.csv" #healthy

    # TO DO - get repertoires for each group from database
    repertoires = [
            ImmuneRepertoireFactory.fromCSVTest(leukemia_path),
            ImmuneRepertoireFactory.fromCSVTest(covid_path),
            ImmuneRepertoireFactory.fromCSVTest(healthy_path),
            ]

    distance_fun = sequenceAligner("BLOSUM62")
    immuneNets = [ simpleBetaDistance(repertoire=repertoire,distance=distance_fun,threshold=0.2) for repertoire in repertoires ]
    
    multiGraphChart(groups,repertoires,immuneNets)
    
