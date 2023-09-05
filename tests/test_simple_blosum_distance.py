import unittest
import pandas as pd
import numpy as np
from src.creation.algoritms.common_methods import *
from src.creation.algoritms.simple_blosum_distance import *
from src.creation.outputStrategies.df_strategy import *
class TestSimpleBlosumDistance(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        path = "test_data/test_clonotypes.csv"
        df = pd.read_csv(path)
        df.name = "some name"
        df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        self.df = df
        self.blosum = BLOSUM62(df_strategy)
        aligner = PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        self.aligner = aligner

    def test_tcr_alig_identical_strings(self):

        result = tcr_alig("GLYYGQ","GLYYGQ",self.aligner)
        self.assertEqual(35, result)  # add assertion here

    def test_tcr_alig_different_strings(self):

        result = tcr_alig("GLAAAQ","GLYYGQ",self.aligner)
        self.assertEqual(11, result)

    def test_blosum_network(self):
        df_net = self.blosum.createGraph(self.df)
        pd.testing.assert_frame_equal(pd.DataFrame(data={'r1': [16,17], 'r2': [15,16]}), df_net)

if __name__ == '__main__':
    unittest.main()
