from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from src.creation.algorithms.common_methods import numerizeTCRSeq
import numpy as np
import math

import pandas as pd
import numpy as np
from src.creation.algorithms.simple_distance import simple_distance
from src.creation.immuneRepertoire import immuneRepertoire
from src.creation.algorithms.common_methods import split_tcr_column
import uuid
from src.creation.utils.pathManager import pathManager

df = pd.read_csv(pathManager().testDataPath / "bigTest.csv").head(1000)
df.name = "aaa"
df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))



class levenshteinDistance:
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
            res = sum([ i != j for i,j in zip(x,ref)])
        else:
            if self.group == False:
                raise ValueError('\"group\" must be set to true when using method in group mode')
            tcr = numerizeTCRSeq(x)
            ceil_all = np.vectorize(math.ceil) 
            res = ceil_all(squareform(pdist(tcr,metric="hamming"))*tcr.shape[1]).astype(float)

        return res

    def __str__(self):
        distance_name = "levenshtein"
        name_string = ""
        if self.group == True:
            name_string += "group_"
        if self.negative == True:
            name_string+= "negative_"
        name_string += distance_name

        return name_string

# if __name__ == "__main__":
#     a = levenshteinDistance()
#     x = "aa"
#     ref = "aaa"
#     print(a.tcr_dist(x,ref))

if __name__ == "__main__":
    df_net = simple_distance(repertoire=immuneRepertoire(df, {uuid.uuid4().hex: len(df) }), distance=levenshteinDistance(group = True))
    print(df_net)

