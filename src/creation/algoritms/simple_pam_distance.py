# return graph (dataframe) based on another df with clonotypes

import numpy as np
import logging
import pandas as pd
from src.creation.algoritms.algorithm import *
from src.creation.algoritms.common_methods import *

class PAM250(algorithm):
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("PAM250")

    def creationAlgorithm(self,clonotypes):
        tcr_npa = clonotypes[["tcra_aa", "tcrb_aa"]].dropna().to_numpy()
        dist_al_trcb = np.zeros(np.shape(tcr_npa)[0] * np.shape(tcr_npa)[0]).reshape(np.shape(tcr_npa)[0],
                                                                                     np.shape(tcr_npa)[0])
        dist_al_trcb += -1
        matrix_len = len(dist_al_trcb)
        for x in range(0, matrix_len):
            logging.info(str(x) + " out of str " + str(matrix_len) + "rows process in triangular similarity matrix")
            for y in range(0 + x, matrix_len):
                dist_al_trcb[x][y] = tcr_alig(tcr_npa[x][0], tcr_npa[y][0],self.aligner) + tcr_alig(tcr_npa[x][1], tcr_npa[y][1],self.aligner)
                dist_al_trcb[y][x] = dist_al_trcb[x][y]


        treshold = 0.8 # todo -> treat as parameter
        for x in range(0, len(dist_al_trcb)):
            self_score = dist_al_trcb[x][x]
            for y in range(0, len(dist_al_trcb)):
                dist_al_trcb[x][y] /= self_score

        dist_al_trcb = np.tril(dist_al_trcb, k=-1)
        matrix_cutoff = np.where(dist_al_trcb > treshold)

        d = {'r1': matrix_cutoff[0], 'r2': matrix_cutoff[1]}
        df_net = pd.DataFrame(data=d)
        df_net.name = clonotypes.name
        return df_net

