import numpy as np
import pandas as pd
from pathlib import Path
import unittest
from src.creation.io_strategies.test_csv_strategy import test_csv_strategy
from src.analysis.methods.publicSimilarity import publicSimilarity


class testPublicSimilarity(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        path = Path(__file__).parent / "test_data/publicTest0.csv"
        path1 = Path(__file__).parent / "test_data/publicTest1.csv"
        path2 = Path(__file__).parent / "test_data/publicTest2.csv"
        pathExpected = Path(__file__).parent / "expected/expected_publicSimilarity"
        self.repertoires = test_csv_strategy().input(path)
        self.repertoires.clones = pd.concat([self.repertoires.clones,test_csv_strategy().input(path1).clones, test_csv_strategy().input(path2).clones],ignore_index=True)
        self.repertoires.clones.name = "testName"
        with open(pathExpected, "r") as f:
            self.expected = eval(f.read())
            self.expected = [ set(i) for i in self.expected]

    def test_publicSimilarity(self):
        result = [ set(i) for i in publicSimilarity(self.repertoires) ]
        self.assertEqual(result, self.expected)
