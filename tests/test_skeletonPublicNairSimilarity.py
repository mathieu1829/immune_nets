import numpy as np
import pandas as pd
from pathlib import Path
import unittest
from src.creation.io_strategies.test_csv_strategy import test_csv_strategy
from analysis.methods.skeletonPublicNairSimilarity import skeleton_similarity


class testSkeletonPublicNairSimilarity(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.path = Path(__file__).parent / "test_data/bigTest.csv"
        self.path2 = Path(__file__).parent / "test_data/bigTest.csv"
        self.repertoires = test_csv_strategy().input(self.path)
        self.repertoires.clones = self.repertoires.clones.append(test_csv_strategy().input(self.path2))

    def test_skeletonPublicNairSimilarity(self):
        skeleton_similarity(self.repertoires)


