import pandas as pd
from src.creation.algorithms.common_methods import *

class csv_strategy():
    def input(path):
        df =  pd.read_csv(path)
        df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        return df
    def output(self, algo, **kwargs):
        immune_net = algo(self, **kwargs)
        immune_net.to_csv(immune_net.df.name+"_net.csv")
        return None 

    
    