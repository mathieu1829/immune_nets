import unittest
from pathlib import Path
from src.factories import ImmuneRepertoireFactory
from src.cluster.comparisonClusterJob import runClusterJob

class TestClusterJob(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        testDataPath = Path(__file__).parent / "test_data"
        covid0Path = testDataPath / "covid_test_clonotypes_0.csv"
        covid1Path = testDataPath / "covid_test_clonotypes_1.csv"
        healthy0Path = testDataPath / "healthy_test_clonotypes_0.csv"
        healthy1Path = testDataPath / "healthy_test_clonotypes_1.csv"
        cls.all_repertoires = {
                                "covid vs healthy": {
                                    "covid": [
                                        ImmuneRepertoireFactory.fromCSV(path=covid0Path, name=f"test covid 0", desc="test repertoire"),
                                        ImmuneRepertoireFactory.fromCSV(path=covid1Path, name=f"test covid 1", desc="test repertoire"),
                                        ],
                                    "healthy": [
                                        ImmuneRepertoireFactory.fromCSV(path=healthy0Path, name=f"test healthy 0", desc="test repertoire"),
                                        ImmuneRepertoireFactory.fromCSV(path=healthy1Path, name=f"test healthy 1", desc="test repertoire"),
                                        ]
                                    },
                                "covid_0 vs covid_1": {
                                    "covid_0": [ImmuneRepertoireFactory.fromCSV(path=covid0Path, name=f"test covid 0", desc="test repertoire")],
                                    "covid_1": [ImmuneRepertoireFactory.fromCSV(path=covid1Path, name=f"test covid 1", desc="test repertoire")]
                                    },
                                "healthy_0 vs healthy_1": {
                                    "healthy_0": [ImmuneRepertoireFactory.fromCSV(path=healthy0Path, name=f"test healthy 0", desc="test repertoire")],
                                    "healthy_1": [ImmuneRepertoireFactory.fromCSV(path=healthy1Path, name=f"test healthy 1", desc="test repertoire")]
                                    }
                                }
    def test_clusterJob(self):
        runClusterJob(self.all_repertoires, 3, True)

if __name__ == '__main__':
    unittest.main()

