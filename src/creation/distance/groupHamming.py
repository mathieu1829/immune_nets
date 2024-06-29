from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from src.creation.algorithms.common_methods import numerizeTCRSeq

import pandas as pd
import numpy as np
from src.creation.algorithms.simple_distance import simple_distance
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.algorithms.common_methods import split_tcr_column
import uuid

df = pd.read_csv("10k_PBMC_5pv2_nextgem_Chromium_Controller_10k_PBMC_5pv2_nextgem_Chromium_Controller_vdj_t_clonotypes.csv").head(1000)
df.name = "aaa"
df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))


class groupHammingDistance:
    def tcr_dist(self,tcr):
        tcr = numerizeTCRSeq(tcr)
        return squareform(pdist(tcr,metric="hamming"))
        

    def __str__(self):
        return "group_hamming"

if __name__ == "__main__":
    df_net = simple_distance(repertoire=immuneRepertoire(df, {uuid.uuid4().hex: len(df) }), distance=groupHammingDistance())
    print(df_net)

