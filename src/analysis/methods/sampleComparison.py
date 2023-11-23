import numpy as np
import pandas as pd

def morisita_horn_similarity(x_stats,y_stats):
    #x/y_stats are placeholders for the time being
    #it will have to be specified whether simpson index should be from tcrb or tcra or both
    #the same with "num of unique"
    #the same with "unique_distribution"
    return float( 2 * (np.array([ x_stats.unique_distribution[i] * y_stats.unique_distribution[i] for i in x_stats.num_of_unique])).sum() ) / float( (x_stats.simpson_index + y_stats.simpson_index) * x_stats.num_of_unique * y_stats.num_of_unique )

def alpha_similarity(stat_list,q):
    stat_list_size = sum([i.unique_distribution.sum() for i in stat_list])
    subunit_abundance = [ ( i.unique_distribution.sum() / stat_list_size)  for i in stat_list] 
    species_abundance = np.array([ [( specimen / species.unique_distribution.sum() ) for specimen in species.unique_distribution] for species in stat_list ])
    return 1 / (np.array([ [ ((subunit_abundance[i]**(q-1)) * species_abundance[i,j]) for j in i.unique_distribution.shape[0] ] for i in range(len(stat_list))]).sum()**(1/(q-1)))
