from pathlib import Path

from immune_nets.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from immune_nets.creation.distance.alignment import sequenceAligner

from immune_nets.mappers import ImmuneNetworkMapper, GraphStatsMapper
from immune_nets.entities import GraphStats
from immune_nets.factories import ImmuneRepertoireFactory

def generateBasicNetworks():
    root_dir = Path(__file__).parent.parent.parent
    test_data_path = root_dir / "tests/test_data"

    leukemia_path_0 = root_dir  / "tests/test_data/leukemia_test_clonotypes_0.csv" # leukemia
    covid_path_0 = root_dir / "tests/test_data/covid_test_clonotypes_0.csv" # covid
    covid_path_1 = root_dir / "tests/test_data/covid_test_clonotypes_1.csv" # covid
    healthy_path_0 = root_dir / "tests/test_data/healthy_test_clonotypes_0.csv" #healthy
    healthy_path_1 = root_dir / "tests/test_data/healthy_test_clonotypes_1.csv" #healthy


    groups = ["healthy", "healthy", "leukemia", "covid", "covid"]
    paths = [healthy_path_0, healthy_path_1, leukemia_path_0, covid_path_0, covid_path_1]
    nums = [0,1,0,0,1]


    distance_fun = sequenceAligner("BLOSUM62")

    statList = []

    for group,path,num in zip(groups,paths,nums):
        repertoire = ImmuneRepertoireFactory.fromCSV(path=path, name=f"{group} repertoire", desc = " ")
         
        net = simpleBetaDistance(repertoire=repertoire,
                                    distance=distance_fun,
                                    threshold=0.3)

        net.name = f"{group} test network_{num}"
        filename = f"{group}_test_network_{num}.pkl"
        filename_csv = f"{group}_test_network_{num}.csv"
        net.graph.to_csv(test_data_path / filename_csv)
        ImmuneNetworkMapper.toPickle(net, test_data_path / filename)
        print(f"Generated: {filename}")
        statList.append(GraphStats(net))



if __name__ == "__main__":
    generateBasicNetworks()
