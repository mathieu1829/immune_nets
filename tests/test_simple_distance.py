import unittest
import pandas as pd
from pathlib import Path
import numpy as np
import uuid
from src.creation.algorithms.common_methods import *
from src.creation.distance.alignment import sequenceAligner
from src.creation.algorithms.simple_distance import *
from src.creation.enums.matrices import *
from src.creation.enums.utils import * 
from src.creation.io_strategies.df_strategy import *
from src.creation.immuneRepertoire import immuneRepertoire
class TestSimpleDistance(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        path = Path(__file__).parent / "test_data/test_clonotypes.csv"
        df = pd.read_csv(path)
        df.name = "MyTest"
        df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        self.df = df
        aligner = PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        self.aligner = aligner

    def test_tcr_alig_identical_strings(self):

        result = sequenceAligner("BLOSUM62").tcr_dist("GLYYGQ","GLYYGQ")
        self.assertEqual(35, result)  # add assertion here

    def test_tcr_alig_different_strings(self):

        result = sequenceAligner("BLOSUM62").tcr_dist("GLAAAQ","GLYYGQ")
        self.assertEqual(11, result)

    def test_simple_distance_networks(self):
        for dist in makeEnumDict(Matrices):
            df_net = simple_distance(repertoire=immuneRepertoire(self.df, {uuid.uuid4().hex: len(self.df) }), distance=sequenceAligner(dist))

            match dist:
                case "PAM250":
                    expected_df = pd.DataFrame(data={'r1': [16,17,17], 'r2': [15,15,16]}) 
                case "PAM30":
                    expected_df = pd.DataFrame(data={'r1': [16], 'r2': [15]}) 
                case _:
                    expected_df = pd.DataFrame(data={'r1': [16,17], 'r2': [15,16]})
            if not df_net.network.equals(expected_df):
                print(f"Error occured while creating matrix: {dist}\n")
                print("Corrupted net:")
                print(df_net.network)
                print("")
            pd.testing.assert_frame_equal(expected_df, df_net.network)

if __name__ == '__main__':
    unittest.main()
