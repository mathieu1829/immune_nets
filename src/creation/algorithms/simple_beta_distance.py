# return graph (dataframe) based on another df with clonotypes

import numpy as np
import logging
import pandas as pd
from src.creation.algorithms.common_methods import *
from src.creation.algorithms.algorithm import *
from src.creation.immuneNetwork import immuneNetwork


@algorithm
def simple_beta_distance(repertoire, distance, threshold = 0.8, **kwargs):
    clonotypes = repertoire.clones
    distanceFun = distance.tcr_dist
    tcr_npa = clonotypes[["tcra_aa", "tcrb_aa"]].dropna().to_numpy()
    dist_al_trcb = np.zeros(np.shape(tcr_npa)[0] * np.shape(tcr_npa)[0]).reshape(np.shape(tcr_npa)[0],
                                                                                 np.shape(tcr_npa)[0])
    if distance.group == False:
        dist_al_trcb += -1
        matrix_len = len(dist_al_trcb)
        for x in range(0, matrix_len):
            logging.info(str(x) + " out of str " + str(matrix_len) + "rows process in triangular similarity matrix")
            for y in range(0 + x, matrix_len):
                dist_al_trcb[x][y] = distanceFun(tcr_npa[x][1], tcr_npa[y][1])
                dist_al_trcb[y][x] = dist_al_trcb[x][y]
    else:
        dist_al_trcb = distanceFun(tcr_npa[:,1])

    # print(dist_al_trcb)

    dist_al_trcb = np.tril(dist_al_trcb, k=-1)
    for i in range(len(dist_al_trcb)):
        for j in range(i,len(dist_al_trcb)):
            dist_al_trcb[i,j] = float('INF')
    matrix_cutoff = np.where(dist_al_trcb < threshold)
    # logging.info(dist_al_trcb)

    d = {'r1': matrix_cutoff[0], 'r2': matrix_cutoff[1]}
    df_net = pd.DataFrame(data=d)
    df_net.name = clonotypes.name
    immuneNet = immuneNetwork(df_net, "simple_distance", np.unique(repertoire.clones["sampleID"].to_numpy()),str(distance) , threshold, len(clonotypes)  ) 

    return immuneNet

