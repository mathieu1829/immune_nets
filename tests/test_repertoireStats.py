import unittest
from pathlib import Path
from immune_nets.entities import RepertoireStats
import pickle 

from immune_nets.creation.algorithms.common_methods import *
from immune_nets.creation.algorithms.simpleDistance import *
from immune_nets.creation.enums.matrices import *
from immune_nets.creation.enums.utils import * 
from immune_nets.factories import ImmuneRepertoireFactory


class TestRepertoireAnalysis(unittest.TestCase):

    
    @classmethod
    def setUpClass(self):
        path = Path(__file__).parent / "test_data/healthy_test_clonotypes_0.csv"
        self.repertoire = ImmuneRepertoireFactory.fromCSVTest(path)

    def listToStr(self, l):
        return [str(i) for i in l]


    def _test_repertoire_analysis(self):
        repertoireStats =  RepertoireStats(self.repertoire)
        repertoireStatsListStr = self.listToStr(repertoireStats.toList())
        # with open("expected_repertoire_stats", "wb") as f:
        #     pickle.dump(repertoireStatsListStr,f)
        with open(Path(__file__).parent / "expected/expected_repertoire_stats", "rb") as f:
            expectedRepertoireStatsListStr = pickle.load(f)
        # print(repertoireStatsListStr)
        # print(expectedRepertoireStatsListStr)
        # print(repertoireStatsListStr == expectedRepertoireStatsListStr)
        # for i,j in zip(repertoireStatsListStr,expectedRepertoireStatsListStr):
        #     print(i == j)
        # print(repertoireStats)
        
        self.assertEqual(repertoireStatsListStr, expectedRepertoireStatsListStr)
    



if __name__ == '__main__':
    unittest.main()
