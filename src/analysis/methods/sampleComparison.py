import numpy as np
import pandas as pd

import uuid
import src.creation.io_strategies.df_strategy 
import src.creation.algorithms.simple_distance 
import src.creation.distance.alignment
from src.creation.algorithms.common_methods import *
from src.creation.distance.alignment import sequenceAligner
from src.creation.algorithms.simple_distance import *
from src.creation.enums.matrices import *
from src.creation.enums.utils import * 
from src.creation.io_strategies.df_strategy import *
from src.creation.immuneRepertoire import immuneRepertoire
from src.analysis.methods.repertoireAnalysis import repertoireAnalysis
df = pd.read_csv("10k_PBMC_5pv2_nextgem_Chromium_Controller_10k_PBMC_5pv2_nextgem_Chromium_Controller_vdj_t_clonotypes.csv").head(1000)
df.name = "aaa"



#formula taken from https://en.wikipedia.org/wiki/Morisita's_overlap_index
def morisita_horn_similarity(x_stats,y_stats):
    #x/y_stats are placeholders for the time being
    #it will have to be specified whether simpson index should be from tcrb or tcra or both
    #the same with "num of unique"
    #the same with "unique_distribution"
    return float( 2 * (np.array([ x_stats.unique_tcrb_distribution[i] * y_stats.unique_tcrb_distribution[i] for i in x_stats.unique_tcrb_distribution.shape[0]])).sum() ) / float( (x_stats.simpson_index_tcrb + y_stats.simpson_index_tcrb) * x_stats.unique_tcrb_distribution.shape[0] * y_stats.unique_tcrb_distribution.shape[0] )

#formula taken from https://en.wikipedia.org/wiki/Alpha_diversity
def alpha_similarity(stat_list,q):
    stat_list_size = sum([i.unique_distribution.sum() for i in stat_list])
    subunit_abundance = [ ( i.unique_distribution.sum() / stat_list_size)  for i in stat_list] 
    species_abundance = np.array([ [( specimen / species.unique_distribution.sum() ) for specimen in species.unique_distribution] for species in stat_list ])
    return 1 / (np.array([ [ ((subunit_abundance[i]**(q-1)) * species_abundance[i,j]) for j in i.unique_distribution.shape[0] ] for i in range(len(stat_list))]).sum()**(1/(q-1)))

if __name__ == "__main__":
    analysis1 = repertoireAnalysis(immuneRepertoire(df, {uuid.uuid4().hex:len(df)}))
    analysis2 = repertoireAnalysis(immuneRepertoire(df, {uuid.uuid4().hex:len(df)}))
