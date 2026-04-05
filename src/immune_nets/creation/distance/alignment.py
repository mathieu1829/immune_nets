from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from immune_nets.creation.algorithms.common_methods import numerizedTCRtoStr
from immune_nets.creation.algorithms.common_methods import numerizeTCRSeq

import pandas as pd
import numpy as np
from immune_nets.creation.algorithms.simpleDistance import simpleDistance
from immune_nets.creation.algorithms.common_methods import split_tcr_column
import uuid
import logging
from immune_nets.creation.utils.pathManager import pathManager

path = pathManager().testDataPath / "healthy_test_clonotypes_0.csv"


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
            x = Seq(x)
            ref = Seq(ref)
            self_score_x = self.__aligner.align(x,x)
            self_score_ref = self.__aligner.align(ref,ref)

            self_score = max(self_score_x.score, self_score_ref.score)
            alignments = self.__aligner.align(x, ref)
            # alignment_min_score = -10 + (min(len(x),len(ref))-1) * -1  
            alignment_min_score = min(0,alignments.score)

            

            res = 1 - ((float(alignments.score)-float(alignment_min_score))/(float(self_score)-float(alignment_min_score))) 
        else:
            if self.group == False:
                raise ValueError('\"group\" must be set to true when using method in group mode')
            tcr = numerizeTCRSeq(x)
            res = squareform(pdist(tcr,metric = lambda u,v : self.__aligner.align(numerizedTCRtoStr(u), numerizedTCRtoStr(v)).score ))
            min_score = res.min()
            max_score = res.max()
            res -= min_score
            res /= (max_score - min_score)


        return (-1 if self.negative == True else 1) * res

    def __str__(self):
        return f"Alignment({self.__matrix})"

# if __name__ == "__main__":
#     a = sequenceAligner(matrix = "BLOSUM62", group = True)
#     b = np.array(["GLYYGQ","GLAAAQ"])
#     print(a.tcr_dist(b))
#     df_net = simpleDistance(repertoire=test_csv_strategy().input(path), distance=sequenceAligner(matrix= "BLOSUM62",group = True))
#     print(df_net)
