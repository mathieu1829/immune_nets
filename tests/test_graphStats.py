import unittest
import pandas as pd
from pathlib import Path
import numpy as np
import uuid
from src.analysis.methods.graphStats import GraphStats
import pickle 

import src.creation.algorithms.simpleDistance 
import src.creation.distance.alignment
from src.creation.algorithms.common_methods import *
from src.creation.distance.alignment import sequenceAligner
from src.creation.algorithms.simpleDistance import *
from src.creation.enums.matrices import *
from src.creation.enums.utils import * 
from src.factories import ImmuneRepertoireFactory
from src.creation.utils.pathManager import pathManager
from src.creation.distance.hamming import hammingDistance


class TestGraphletComposition(unittest.TestCase):

    
    @classmethod
    def setUpClass(self):
        self.path = Path(__file__).parent / "test_data/healthy_test_clonotypes_0.csv"
        self.df_net = simpleDistance(repertoire=ImmuneRepertoireFactory.fromCSVTest(self.path), distance = hammingDistance(group=True))

    def listToStr(self, l):
        return [str(i) for i in l]


    def test_graphletComposition(self):
        graphStats = GraphStats(self.df_net)
        graphStatsList = self.listToStr(graphStats.toList())
        # with open("expected_graphStats","wb") as f:
        #     pickle.dump(graphStatsList,f)
        with open(Path(__file__).parent / "expected/expected_graphStats", "rb") as f:
            expectedGraphStatList = pickle.load(f)
        # for i,j in zip(graphStatsList, expectedGraphletComposition):
        #     print(i == j)

        # print(graphStatsList[6])
        # print(expectedGraphletComposition[6])

        
        self.assertEqual(graphStatsList, expectedGraphStatList)
    



if __name__ == '__main__':
    unittest.main()
