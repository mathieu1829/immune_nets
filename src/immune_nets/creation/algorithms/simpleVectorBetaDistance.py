# return graph (dataframe) based on another df with clonotypes

import numpy as np
from immune_nets.creation.algorithms.common_methods import *
from immune_nets.creation.algorithms.algorithm import *
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from immune_nets.entities import ImmuneNetwork

from immune_nets.models import Network

@algorithm
def simpleVectorBetaDistance( repertoire, distance, threshold = None, **kwargs):
    clonotypes = repertoire.clones
    distanceFun = distance.tcr_dist

    if distance.group == True:
        raise ValueError("Vector distance algorithm can't handle \"group\" distance functions")


    b_tcr = clonotypes[["tcrb_aa"]].dropna();
    tcr_npa = clonotypes["tcrb_aa"].dropna().to_numpy()
    dist_al_trcb = np.zeros(np.shape(tcr_npa)[0] * np.shape(tcr_npa)[0]).reshape(np.shape(tcr_npa)[0],
                                                                                 np.shape(tcr_npa)[0])
    dist_al_trcb += -1

    unique_amino_acids = np.array(['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V',''])

    beta_profile = []

    for seq in range(tcr_npa.shape[0]):
        beta_seq = list(tcr_npa[seq])
        for idx,amino in enumerate(beta_seq):
            if idx == len(beta_profile):
                beta_profile.append(np.zeros(unique_amino_acids.shape[0]))
                beta_profile[-1][-1] += seq
            beta_profile[idx][np.where(unique_amino_acids == amino)]+=1
        for i in range(len(beta_profile) - (len(beta_seq))):
            beta_profile[len(beta_seq) - 1 + i][-1]+=1



    consensus_beta_arr = np.array([unique_amino_acids[beta_profile[peptide].argmax()] for peptide in range(len(beta_profile))])

    consensus_beta_seq = ""
    for i in range(len(beta_profile)):
        consensus_beta_seq+=consensus_beta_arr[i]
    b_tcr['beta_closest_to_consensus'] = b_tcr['tcrb_aa'].apply(lambda x: distanceFun(x,consensus_beta_seq))


    closest_beta = b_tcr[['tcrb_aa','beta_closest_to_consensus']].sort_values('beta_closest_to_consensus').dropna().to_numpy()[:6,0]

    for idx,amino in enumerate(closest_beta):
        b_tcr[f'b{idx}'] = b_tcr['tcrb_aa'].apply(lambda x : distanceFun(x,amino))

    dist_mat = pd.DataFrame(
        squareform(pdist(b_tcr[["b0", "b1", "b2", "b3", "b4", "b5"]])),
        columns = b_tcr.index,
        index = b_tcr.index
    ) 
    
    threshold = np.nanmean(dist_mat.to_numpy()) / 4 if threshold is None else threshold
    # threshold = np.mean(dist_mat.fillna(0).to_numpy()) / 4 if threshold is None else threshold #produces smaller threshold
    dist_mat = np.tril(dist_mat, k=-1)
    dist_mat = np.where(dist_mat == 0, np.inf, dist_mat)
    # threshold = np.nanmean(dist_mat) / 4 if threshold is None else threshold #produces even smaller threshold
    # print(dist_mat)
    # print(f"threshold: {threshold}")
    matrix_cutoff = np.where(dist_mat < threshold)

    d = {'r1': matrix_cutoff[0], 'r2': matrix_cutoff[1]}
    df_net = pd.DataFrame(data=d)
    immuneNet = ImmuneNetwork(graph=df_net,
                              method="simpleVectorBetaDistance",
                              sampleId=repertoire.repertoire_id,
                              distanceFun=str(distance),
                              threshold=threshold, 
                              sampleSize=len(clonotypes),
                              proportions=clonotypes["proportion"].to_numpy()
                              )
    return immuneNet


