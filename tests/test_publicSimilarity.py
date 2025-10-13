import numpy as np
import pandas as pd
from pathlib import Path
import unittest
from src.analysis.methods.publicSimilarity import publicSimilarity
from src.creation.immuneRepertoire import ImmuneRepertoire


class testPublicSimilarity(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        path = Path(__file__).parent / "test_data/publicTest0.csv"
        path1 = Path(__file__).parent / "test_data/publicTest1.csv"
        path2 = Path(__file__).parent / "test_data/publicTest2.csv"
        pathExpected = Path(__file__).parent / "expected/expected_publicSimilarity"
        self.repertoires = ImmuneRepertoire.fromCSVTest(path)
        self.repertoires.clones = pd.concat([self.repertoires.clones,ImmuneRepertoire.fromCSVTest(path1).clones, ImmuneRepertoire.fromCSVTest(path2).clones],ignore_index=True)
        self.repertoires.clones.name = "testName"
        with open(pathExpected, "r") as f:
            self.expected = eval(f.read())
            self.expected = [ set(i) for i in self.expected]

    def _test_publicSimilarity(self):
        result = [ set(i) for i in publicSimilarity(self.repertoires) ]
        self.assertEqual(result, self.expected)
