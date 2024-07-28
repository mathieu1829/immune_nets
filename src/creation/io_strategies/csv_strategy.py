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
        sampleIDs = {}
        dfSize = 0;
        names = []
        for path in paths:
            sample_df =  pd.read_csv(path)
            df_list.append(sample_df)
            newUUID = path.split('/')[-1][:32] 
            dfSize+=len(sample_df)
            names.append(path.split('/')[-1])
            if checkUUID(newUUID):
                sampleIDs[newUUID] = dfSize
                # df_list[-1]["sampleID"] = np.array([newUUID] * df_list[-1].shape[0])
            else:
                newUUID = uuid.uuid4().hex
                sampleIDs[newUUID] = dfSize
                newPath = path.split("/")
                newPath[-1] =  f'{newUUID}_{newPath[-1]}'
                newPath = "/".join(newPath)
                os.rename(path,newPath)
            df_list[-1]["sampleID"] = np.array([newUUID] * df_list[-1].index)
        df = pd.concat(df_list)
        df.name = "_".join(names)
        df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        return immuneRepertoire(clones = df, sampleIDs = sampleIDs)
    def output(self, algo, **kwargs):
        immune_net = algo(**kwargs)
        immune_net.network.to_csv(immune_net.network.name+"_net.csv")
        with open(f"{immune_net.network.name}_sampleIDs") as f:
            f.write(repr(immune_net.sampleIDs))
        return None 

    
    
