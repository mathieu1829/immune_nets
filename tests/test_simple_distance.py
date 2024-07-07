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
from src.creation.io_strategies.test_csv_strategy import *
from src.creation.immuneRepertoire import immuneRepertoire
class TestSimpleDistance(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.path = Path(__file__).parent / "test_data/test_clonotypes.csv"
        

    def test_tcr_alig_identical_strings(self):

        result = sequenceAligner("BLOSUM62").tcr_dist("GLYYGQ","GLYYGQ")
        # modelValue = 35
        modelValue = 0.0
        if modelValue != result:
            print(f"model value was {modelValue} but the actual result was {result}")
        self.assertEqual(modelValue, result)  # add assertion here

    def test_tcr_alig_different_strings(self):
        result = sequenceAligner("BLOSUM62").tcr_dist("GLAAAQ","GLYYGQ")
        # modelValue = 11
        modelValue = 0.6857142857142857
        if modelValue != result:
            print(f"model value was {modelValue} but the actual result was {result}")
        self.assertEqual(modelValue, result)

    def test_simple_distance_networks(self):
        for dist in makeEnumDict(Matrices):
            print(f"testing distance: {dist}")
            df_net = simple_distance(repertoire=test_csv_strategy().input(self.path), distance=sequenceAligner(dist))

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
            print(f"test for {dist} completed succesfully")

if __name__ == '__main__':
    unittest.main()
