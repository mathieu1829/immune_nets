import unittest
import pandas as pd
from pathlib import Path
import numpy as np
import uuid
from src.analysis.methods.repertoireAnalysis import repertoireAnalysis
import pickle 

import src.creation.algorithms.simpleDistance 
import src.creation.distance.alignment
from src.creation.algorithms.common_methods import *
from src.creation.distance.alignment import sequenceAligner
from src.creation.algorithms.simpleDistance import *
from src.creation.enums.matrices import *
from src.creation.enums.utils import * 
from src.creation.immuneRepertoire import ImmuneRepertoire
from src.creation.utils.pathManager import pathManager
from src.creation.distance.hamming import hammingDistance


class TestRepertoireAnalysis(unittest.TestCase):

    
    @classmethod
    def setUpClass(self):
        path = Path(__file__).parent / "test_data/healthy_test_clonotypes_0.csv"
        self.repertoire = ImmuneRepertoire.fromCSVTest(path)

    def listToStr(self, l):
        return [str(i) for i in l]


    def _test_repertoire_analysis(self):
        repertoireStats =  repertoireAnalysis(self.repertoire)
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
