import numpy as np
import pandas as pd
from pathlib import Path
import unittest
from src.creation.io_strategies.test_csv_strategy import test_csv_strategy
from src.analysis.methods.skeletonPublicNairSimilarity import skeletonPublicNairSimilarity


class testSkeletonPublicNairSimilarity(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        path = Path(__file__).parent / "test_data/publicTest0.csv"
        path1 = Path(__file__).parent / "test_data/publicTest1.csv"
        path2 = Path(__file__).parent / "test_data/publicTest2.csv"
        pathExpected = Path(__file__).parent / "expected/expected_skeletonPublicNairSimilarity"
        self.repertoires = test_csv_strategy().input(path)
        self.repertoires.clones = pd.concat([self.repertoires.clones,test_csv_strategy().input(path1).clones, test_csv_strategy().input(path2).clones],ignore_index=True)
        self.repertoires.clones.name = "testName"

        with open(pathExpected, "r") as f:
            self.expected = eval(f.read())
            self.expected = [ set(i) for i in self.expected]
        

    def test_skeletonPublicNairSimilarity(self):

        # tcrb_lenghts = test_csv_strategy().input(self.path).clones["tcrb_aa"].dropna().apply(lambda x : len(x))
        # tcra_lenghts = test_csv_strategy().input(self.path).clones["tcra_aa"].dropna().apply(lambda x : len(x))
        # print(f"""
        # tcrb:
        #     mean lenght: {np.mean(tcrb_lenghts)}
        #     std: {np.std(tcrb_lenghts)}
        #     max: {np.max(tcrb_lenghts)}
        #     min: {np.min(tcrb_lenghts)}
        # """
        # )
        # print(f"""
        # tcra:
        #     mean lenght: {np.mean(tcra_lenghts)}
        #     std: {np.std(tcra_lenghts)}
        #     max: {np.max(tcra_lenghts)}
        #     min: {np.min(tcra_lenghts)}
        # """
        # )
        result = [ set(i) for i in skeletonPublicNairSimilarity(self.repertoires) ]
        self.assertEqual(result, self.expected)



