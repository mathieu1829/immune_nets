import pandas as pd
import numpy as np
import uuid
import os
from src.creation.algorithms.common_methods import *
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.utils.checkUUID import checkUUID

class csv_strategy():
    def input(self,paths):
        if type(paths) == str:
            paths = [paths]
        df_list = []
        sampleIDs = []
        for path in paths:
            sample_df =  pd.read_csv(path)
            df_list.append(sample_df)
            hardcodeID = path.split('/')[-1][:32] 
            if checkUUID(hardcodeID):
                sampleIDs.append(hardcodeID)
            else:
                newUUID = uuid.uuid4().hex
                sampleIDs.append(newUUID)
                newPath = path.split("/")
                newPath[-1] =  f'{newUUID}_{newPath[-1]}'
                newPath = "/".join(newPath)
                os.rename(path,newPath)
        df = pd.concat(df_list)
        df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        return immuneRepertoire(clones = df, sampleIDs = sampleIDs)
    def output(self, algo, **kwargs):
        immune_net = algo(**kwargs)
        immune_net.network.to_csv(immune_net.network.name+"_net.csv")
        with open(f"{immune_net.network.name}_sampleIDs") as f:
            f.write(repr(immune_net.sampleIDs))
        return None 

    
    
