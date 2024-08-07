import unittest
import pandas as pd
from pathlib import Path
import numpy as np
import uuid
from src.analysis.methods.graphlet_composition import graphletComposition

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

    def test_graphlet_coposition(self):
        with open(Path(__file__).parent / "expected/expected_graphlet_composition","r") as f:
            expected_graphlet_composition = f.read()
        graph_stats = graphletComposition(self.df_net)
        print(graph_stats.vertice_num)
        print(graph_stats.isolated_vertices)
        self.assertEqual(repr(graph_stats.toList()),expected_graphlet_composition)
    



if __name__ == '__main__':
    unittest.main()
