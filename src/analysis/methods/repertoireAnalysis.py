import numpy as np
import pandas as pd
from src.creation.io_strategies.csv_strategy import csv_strategy


df = pd.read_csv("10k_PBMC_5pv2_nextgem_Chromium_Controller_10k_PBMC_5pv2_nextgem_Chromium_Controller_vdj_t_clonotypes.csv")

def repertoireAnalysis(repertoire):
    tcra = np.unique(repertoire['tcra_aa'].dropna().to_numpy())
    tcrb = np.unique(repertoire['tcrb_aa'].dropna().to_numpy())
    all_tcr = np.unique(np.concatenate([tcra,tcrb]))

    num_of_tcra = tcra.shape[0]
    num_of_tcrb = tcrb.shape[0]
    num_of_all_tcr = all_tcr.shape[0]

    unique_tcra_distribution = np.array([ (tcra == seq).sum()  for seq in tcra ])
    unique_tcrb_distribution = np.array([ (tcrb == seq).sum()  for seq in tcrb ])
    unique_all_tcr_distribution = np.array([ (all_tcr == seq).sum()  for seq in all_tcr ])
    
    simpson_index_tcra = (unique_tcra_distribution**2).sum() 
    simpson_index_tcrb = (unique_tcrb_distribution**2).sum() 
    simpson_index_all_tcr = (unique_all_tcr_distribution**2).sum() 
    transform = lambda x : x * np.log(x)
    shannon_index_tcra = transform(unique_tcra_distribution).sum()
    shannon_index_tcrb = transform(unique_tcrb_distribution).sum()
    shannon_index_all_tcr = transform(unique_all_tcr_distribution).sum()
    
    print(simpson_index_tcrb)

if __name__ == "__main__":
    repertoireAnalysis(csv_strategy().input("10k_PBMC_5pv2_nextgem_Chromium_Controller_10k_PBMC_5pv2_nextgem_Chromium_Controller_vdj_t_clonotypes.csv"))
