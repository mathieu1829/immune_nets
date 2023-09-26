import pandas as pd
from src.creation.algorithms.common_methods import *

def csv_strategy(path):
    df =  pd.read_csv(path)
    df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
    df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
    return df