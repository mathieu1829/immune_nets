from scipy.spatial.distance import hamming
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from src.creation.algorithms.common_methods import numerizeTCRSeq
import numpy as np

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


class hammingDistance:
    def __init__(self, group = False, negative = False):
        self.group = group
        self.negative = negative

    def tcr_dist(self,x,ref=None):
        if x is None :
            raise ValueError('\"x\" parameter must not be none')

        if ref is not None:
            if self.group == True:
                raise ValueError('\"group\" must not be True, when using method in single comparison mode')
            maxLen = max(len(x), len(ref))
            x = x.ljust(maxLen)
            ref = ref.ljust(maxLen)
            res = hamming(list(x), list(ref))
        else:
            if self.group == False:
                raise ValueError('\"group\" must be set to true when using method in group mode')
            tcr = numerizeTCRSeq(x)
            res = squareform(pdist(tcr,metric="hamming"))

        return res
        

    def __str__(self):
        distance_name = "hamming"
        name_string = ""
        if self.group == True:
            name_string += "group_"
        if self.negative == True:
            name_string+= "negative_"
        name_string += distance_name

        return name_string

if __name__ == "__main__":
    df_net = simple_distance(repertoire=immuneRepertoire(df, {uuid.uuid4().hex: len(df) }), distance=hammingDistance(group = False))
    print(df_net)
