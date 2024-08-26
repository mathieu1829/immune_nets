import unittest
import pandas as pd
from pathlib import Path
import numpy as np
import uuid
from src.analysis.methods.graphlet_composition import graphletComposition
import pickle 

import src.creation.algorithms.simple_distance 
import src.creation.distance.alignment
from src.creation.algorithms.common_methods import *
from src.creation.distance.alignment import sequenceAligner
from src.creation.algorithms.simple_distance import *
from src.creation.enums.matrices import *
from src.creation.enums.utils import * 
from src.creation.io_strategies.test_csv_strategy import *
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.utils.pathManager import pathManager
from src.creation.distance.hamming import hammingDistance


class TestGraphletComposition(unittest.TestCase):

    
    @classmethod
    def setUpClass(self):
        self.path = Path(__file__).parent / "test_data/bigTest.csv"
        self.df_net = simple_distance(repertoire=test_csv_strategy().input(self.path), distance = hammingDistance(group=True))

    def listToStr(self, l):
        return [str(i) for i in l]


    def test_graphlet_composition(self):
        graphStats = graphletComposition(self.df_net)
        graphStatsList = self.listToStr(graphStats.toList())
        # with open("expected_graphlet_composition","wb") as f:
        #     pickle.dump(graphStatsList,f)
        with open(Path(__file__).parent / "expected/expected_graphlet_composition", "rb") as f:
            expectedGraphletComposition = pickle.load(f)
        # for i,j in zip(graphStatsList, expectedGraphletComposition):
        #     print(i == j)

        # print(graphStatsList[6])
        # print(expectedGraphletComposition[6])

        
        self.assertEqual(graphStatsList, expectedGraphletComposition)
    



if __name__ == '__main__':
    unittest.main()
