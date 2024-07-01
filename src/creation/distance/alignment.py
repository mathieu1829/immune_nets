from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from src.creation.algorithms.common_methods import numerizedTCRtoStr
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

class sequenceAligner:
    def setMatrix(self,matrix):
        self.__matrix = matrix
        self.__aligner = PairwiseAligner()
        self.__aligner.substitution_matrix = substitution_matrices.load(matrix)
        self.__aligner.open_gap_score = -10;
        self.__aligner.extend_gap_score = -1;

    def getMatrix(self):
        return self.__matrix

    def __init__(self, matrix,switch_self_score=False, group = False, negative = False):
        self.setMatrix(matrix)
        self.switch_self_score = switch_self_score
        self.group = group
        self.negative = negative

        
        

    def tcr_dist(self,x, ref = None):
        if x is None :
            raise ValueError('\"x\" parameter must not be none')

        if ref is not None:
            if self.group == True:
                raise ValueError('\"group\" must not be True, when using method in single comparison mode')
            self_score_x = self.__aligner.align(x,x)
            self_score_ref = self.__aligner.align(ref,ref)

            self_score = max(self_score_x.score, self_score_ref.score)
            alignments = self.__aligner.align(x, ref)

            res = 1 - (float(alignments.score) / float(self_score)) 
        else:
            if self.group == False:
                raise ValueError('\"group\" must be set to true when using method in group mode')
            tcr = numerizeTCRSeq(x)
            res = squareform(pdist(tcr,metric = lambda u,v : self.__aligner.align(numerizedTCRtoStr(u), numerizedTCRtoStr(v)).score ))
            min_score = res.min()
            max_score = res.max()
            res -= min_score
            res /= (max_score - min_score)


        return res

    def __str__(self):
        return self.__matrix

if __name__ == "__main__":
    a = sequenceAligner(matrix = "BLOSUM62", group = True)
    b = np.array(["GLYYGQ","GLAAAQ"])
    print(a.tcr_dist(b))
    df_net = simple_distance(repertoire=immuneRepertoire(df, {uuid.uuid4().hex: len(df) }), distance=sequenceAligner(matrix= "BLOSUM62",group = True))
    print(df_net)
