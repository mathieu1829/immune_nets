import pandas as pd
import numpy as np
import uuid
import os
from src.creation.algorithms.common_methods import *
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.utils.checkUUID import checkUUID

class test_csv_strategy():
    def input(self,path):

        df =  pd.read_csv(path)
        df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        
        return immuneRepertoire(clones = df, sampleIDs = uuid.uuid4().hex)
    def output(self, algo, **kwargs):
        immune_net = algo(**kwargs)
        immune_net.network.to_csv(immune_net.network.name+"_net.csv")
        with open(f"{immune_net.network.name}_sampleIDs") as f:
            f.write(repr(immune_net.sampleIDs))
        return None 

    
    
