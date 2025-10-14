import numpy as np
import pandas as pd
from pathlib import Path
import unittest
from src.analysis.methods.privateSimilarity import privateSimilarity
from src.factories import ImmuneRepertoireFactory


class testPrivateSimilarity(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        path = Path(__file__).parent / "test_data/publicTest0.csv"
        path1 = Path(__file__).parent / "test_data/publicTest1.csv"
        path2 = Path(__file__).parent / "test_data/publicTest2.csv"
        pathExpected = Path(__file__).parent / "expected/expected_privateSimilarity"
        self.repertoires = ImmuneRepertoireFactory.fromCSVTest(path)
        self.repertoires.clones = pd.concat([self.repertoires.clones,ImmuneRepertoireFactory.fromCSVTest(path1).clones, ImmuneRepertoireFactory.fromCSVTest(path2).clones],ignore_index=True)
        self.repertoires.clones.name = "testName"

        with open(pathExpected, "r") as f:
            self.expected = eval(f.read())
            # self.expected = [ [ set(cluster) for cluster in sample ] for sample in self.expected]
            self.expected = set([ i for sample in self.expected for cluster in sample for i in cluster])
        

    def _test_privateSimilarity(self):

        result = set([ i for sample in privateSimilarity(self.repertoires) for cluster in sample for i in cluster]) 
        self.assertEqual(result, self.expected)
