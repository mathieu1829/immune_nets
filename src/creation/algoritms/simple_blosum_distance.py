# return graph (dataframe) based on another df with clonotypes

import numpy as np
from Bio import pairwise2
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import logging
import pandas as pd

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")


def tcr_alig(x, ref):
    if x != None:
        alignments = pairwise2.align.localds(x, ref, aligner.substitution_matrix, -10, -1)
        return alignments[0].score


def create_BLOSUM62_network(df):
    tcr_npa = df[["tcra_aa", "tcrb_aa"]].dropna().to_numpy()
    dist_al_trcb = np.zeros(np.shape(tcr_npa)[0] * np.shape(tcr_npa)[0]).reshape(np.shape(tcr_npa)[0],
                                                                                 np.shape(tcr_npa)[0])
    dist_al_trcb += -1
    matrix_len = len(dist_al_trcb)
    for x in range(0, matrix_len):
        logging.info(str(x) + " out of str " + str(matrix_len) + "rows process in triangular similarity matrix")
        for y in range(0 + x, matrix_len):
            dist_al_trcb[x][y] = tcr_alig(tcr_npa[x][0], tcr_npa[y][0]) + tcr_alig(tcr_npa[x][1], tcr_npa[y][1])
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
    return df_net

