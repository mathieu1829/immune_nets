from pathlib import Path

from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.distance.alignment import sequenceAligner

from src.mappers import ImmuneNetworkMapper, GraphStatsMapper
from src.analysis.methods.graphStats import GraphStats
from src.factories import ImmuneRepertoireFactory

def generateBasicNetworks():
    root_dir = Path(__file__).parent.parent.parent
    test_data_path = root_dir / "tests/test_data"

    leukemia_path = root_dir  / "tests/test_data/leukemia_test_clonotypes.csv" # leukemia
    covid_path = root_dir / "tests/test_data/covid_test_clonotypes.csv" # covid
    healthy_path = root_dir / "tests/test_data/healthy_test_clonotypes_1.csv" #healthy


    groups = ["healthy", "leukemia", "covid"]
    paths = [healthy_path, leukemia_path, covid_path]


    distance_fun = sequenceAligner("BLOSUM62")

    statList = []

    for group,path in zip(groups,paths):
        repertoire = ImmuneRepertoireFactory.fromCSV(path=path, name=f"{group} repertoire", desc = " ")
         
        net = simpleBetaDistance(repertoire=repertoire,
                                    distance=distance_fun,
                                    threshold=0.3)

        net.name = f"{group} test network"
        filename = f"{group}_test_network.pkl"
        filename_csv = f"{group}_test_network.csv"
        net.graph.to_csv(test_data_path / filename_csv)
        ImmuneNetworkMapper.toPickle(net, test_data_path / filename)
        print(f"Generated: {filename}")
        statList.append(GraphStats(net))



if __name__ == "__main__":
    generateBasicNetworks()
