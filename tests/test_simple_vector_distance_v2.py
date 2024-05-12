import unittest
import pandas as pd
from pathlib import Path
import numpy as np
import os
import uuid
from src.creation.algorithms.common_methods import *
from src.creation.algorithms.simple_vector_distance_v2 import *
from src.creation.enums.matrices import *
from src.creation.enums.utils import * 
from src.creation.distance.alignment import sequenceAligner
from src.creation.io_strategies.df_strategy import *
class TestSimpleVectorDistance(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        path = Path(__file__).parent / "test_data/test_clonotypes.csv"
        df = pd.read_csv(path)
        df.name = "MyTest"
        df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        self.df = df


    def test_simple_vector_network(self):
        for dist in makeEnumDict(Matrices):
            df_net = simple_vector_distance_v2(repertoire=immuneRepertoire(self.df, {uuid.uuid4().hex: len(self.df) }), distance=sequenceAligner(dist))
            expected_df = pd.DataFrame(data={'r1': [16,17], 'r2': [15,16]})
            match dist:
                case "PAM250":
                    expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_PAM250_vector_v2.csv")[["r1","r2"]]
                case "PAM30":
                    expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_PAM30_vector_v2.csv")[["r1","r2"]]
                case "PAM70":
                    expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_PAM70_vector_v2.csv")[["r1","r2"]]
                case "BLOSUM45":
                    expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_BLOSUM45_vector_v2.csv")[["r1","r2"]]
                case "BLOSUM50":
                    expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_BLOSUM50_vector_v2.csv")[["r1","r2"]]
                case "BLOSUM62":
                    expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_BLOSUM62_vector_v2.csv")[["r1","r2"]]
                case "BLOSUM80":
                    expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_BLOSUM80_vector_v2.csv")[["r1","r2"]]
                case "BLOSUM90":
                    expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_BLOSUM90_vector_v2.csv")[["r1","r2"]]
                case _:
                    expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_PAM250_vector.csv")[["r1","r2"]]
            if not df_net.network.equals(expected_df) and not (expected_df.empty and df_net.network.empty):
                print(f"Error occured while creating matrix: {dist}\n")
                print("Corrupted net:")
                print(df_net.network)
                print("Good net:")
                print(expected_df)
                # new_expected_df = df_net
                # found = False
                # for file in os.listdir("tests/expected") :
                #     test_df = pd.read_csv(Path(__file__).parent / "expected" / file )[["r1","r2"]]
                #     if df_net.equals(test_df):
                #         new_expected_df = test_df
                #         print("new expected_df exists: "+file)
                #         found = True
                #         break
                # if not found:        
                #     print("new expected_df does not exits, output to file")
                #     new_expected_df.to_csv(f"expected_{dist}_vector_v2.csv")
                print("")
                # df_net.network.to_csv(f"expected_{dist}_vector_v2.csv")
            if not (expected_df.empty and df_net.network.empty):
                pd.testing.assert_frame_equal(expected_df, df_net.network)

if __name__ == '__main__':
    unittest.main()
