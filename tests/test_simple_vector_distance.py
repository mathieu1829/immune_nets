import unittest
import pandas as pd
from pathlib import Path
import numpy as np
from src.creation.algorithms.common_methods import *
from src.creation.algorithms.simple_vector_distance import *
from src.creation.enums.matrices import *
from src.creation.enums.utils import * 
from src.creation.outputStrategies.df_strategy import *
class TestSimpleVectorDistance(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        path = Path(__file__).parent / "test_data/test_clonotypes.csv"
        df = pd.read_csv(path)
        df.name = "MyTest"
        df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        self.df = df
        self.simple_vector = simple_vector_distance(df_strategy)


    def test_simple_vector_network(self):
        for dist in ["BLOSUM62"]: #makeEnumDict(Matrices):
            df_net = self.simple_vector.createGraph(clonotypes=self.df, matrix=getattr(Matrices,dist))
            expected_df = pd.DataFrame(data={'r1': [16,17], 'r2': [15,16]})
            # match dist:
            #     case "PAM250":
            #         expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_PAM250_vector.csv")[["r1","r2"]]
            #     case "PAM30":
            #         expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_PAM30_vector.csv")[["r1","r2"]]
            #     case "PAM70":
            #         expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_PAM70_vector.csv")[["r1","r2"]]
            #     case "BLOSUM45":
            #         expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_BLOSUM45_vector.csv")[["r1","r2"]]
            #     case "BLOSUM50":
            #         expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_BLOSUM50_vector.csv")[["r1","r2"]]
            #     case "BLOSUM62":
            #         expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_BLOSUM62_vector.csv")[["r1","r2"]]
            #     case "BLOSUM80":
            #         expected_df = pd.read_csv(Path(__file__).parent / "expected/expected_BLOSUM80_vector.csv")[["r1","r2"]]
            #     case _:
            #         expected_df = pd.DataFrame(data={'r1': [16,17], 'r2': [15,16]})
            if not df_net.equals(expected_df):
                print(f"Error occured while creating matrix: {dist}\n")
                print("Corrupted net:")
                print(df_net)
                print("Good net:")
                print(expected_df)
                df_net.to_csv(f"expected_{dist}_vector.csv")
                print("")
            pd.testing.assert_frame_equal(expected_df, df_net)

if __name__ == '__main__':
    unittest.main()