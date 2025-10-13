# return graph (dataframe) based on another df with clonotypes

import numpy as np
import logging
import pandas as pd
from src.creation.algorithms.common_methods import *
from src.creation.algorithms.algorithm import *
from src.creation.immuneNetwork import immuneNetwork


@algorithm
def simpleAlphaDistance(repertoire, distance, threshold = 0.8, **kwargs):
    clonotypes = repertoire.clones
    distanceFun = distance.tcr_dist
    tcr_npa = clonotypes[["tcra_aa", "tcrb_aa"]].dropna().to_numpy()
    dist_al_trcb = np.full((np.shape(tcr_npa)[0], np.shape(tcr_npa)[0]),np.inf)

    if distance.group == False:
        matrix_len = len(dist_al_trcb)
        for y in range(0, matrix_len):
            logging.info(str(y) + " out of str " + str(matrix_len) + "rows process in triangular similarity matrix")
            for x in range(0, y+1):
                if x != y:
                    dist_al_trcb[y][x] = distanceFun(tcr_npa[x][0], tcr_npa[y][0]) 
    else:
        dist_al_trcb = distanceFun(tcr_npa[:,0])

    matrix_cutoff = np.where(dist_al_trcb < threshold)

    d = {'r1': matrix_cutoff[0], 'r2': matrix_cutoff[1]}
    df_net = pd.DataFrame(data=d)
    df_net.name = clonotypes.name
    immuneNet = immuneNetwork(graph=df_net,
                              method="simpleAlphaDistance",
                              sampleId=repertoire.repertoire_id,
                              distanceFun=str(distance),
                              threshold=threshold, 
                              sampleSize=len(clonotypes)
                              )

    return immuneNet

