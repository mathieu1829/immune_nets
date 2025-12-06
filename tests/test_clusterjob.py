import unittest
from pathlib import Path
from src.factories import ImmuneRepertoireFactory
from src.cluster.clusterJob import runClusterJob

class TestClusterJob(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        testDataPath = Path(__file__).parent / "test_data"
        covid0Path = testDataPath / "covid_test_clonotypes_0.csv"
        covid1Path = testDataPath / "covid_test_clonotypes_1.csv"
        cls.all_repertoires = {"covid0 vs covid1": {
                                    "covid0": [ImmuneRepertoireFactory.fromCSV(path=covid0Path, name=f"test covid 0", desc="test repertoire")],
                                    "covid1": [ImmuneRepertoireFactory.fromCSV(path=covid1Path, name=f"test covid 1", desc="test repertoire")]
                                    }
                                }
    def test_clusterJob(self):
        runClusterJob(self.all_repertoires, 3, True)

if __name__ == '__main__':
    unittest.main()

