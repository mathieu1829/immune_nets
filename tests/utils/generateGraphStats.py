from pathlib import Path

from src.factories import ImmuneNetworkFactory
from src.mappers import GraphStatsMapper
from src.analysis.methods.graphStats import GraphStats

def generateGraphStats():
    root_dir = Path(__file__).parent.parent.parent
    test_data_path = root_dir / "tests/test_data"

    leukemia_network_path = test_data_path  / "leukemia_test_network.pkl" # leukemia
    covid_network_path = test_data_path / "covid_test_network.pkl" # covid
    healthy_network_path = test_data_path / "healthy_test_network.pkl" #healthy

    networkPaths = [healthy_network_path, leukemia_network_path, covid_network_path]
    
    statList = [ ImmuneNetworkFactory.fromPickle(path) for path in networkPaths]

    statList = [ GraphStats(net) for net in statList]

    


    GraphStatsMapper.toDataframe(statList).to_pickle(test_data_path / "graphStats.pkl")

    import pandas as pd

    print("pickle:")
    df = pd.read_pickle(test_data_path / "graphStats.pkl")
    print(df)

if __name__ == "__main__":
    generateGraphStats()
