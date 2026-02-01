import numpy as np
from scipy.stats import wasserstein_distance
from src.analysis.statDistanceTypes import PairwiseStatDistance

from src.analysis.distributionInference import DistributionIdentifier, NaiveIdentifier
from src.entities import GraphStats

class PairwiseDistributionDistance(PairwiseStatDistance):
    def __init__(self, distributionName, identifier: DistributionIdentifier = NaiveIdentifier()):
        self.distributionName = distributionName
        self.identifier = identifier

    def pmf_mean(self, dist):
        values = np.array(list(dist.keys()))
        probs  = np.array(list(dist.values()))
        return np.sum(values * probs)

    def pmf_min(self, dist):
        values = np.array(list(dist.keys()))
        return np.min(values)

    def pmf_max(self, dist):
        values = np.array(list(dist.keys()))
        return np.max(values)


    def pmf_var(self, dist):
        values = np.array(list(dist.keys()))
        probs  = np.array(list(dist.values()))
        mu = self.pmf_mean(dist)
        return np.sum(probs * (values - mu)**2)

    def pmf_std(self, dist):
        return np.sqrt(self.pmf_var(dist)) 

    def pmf_gini(self, dist):
        values = np.array(list(dist.keys()))
        probs  = np.array(list(dist.values()))
        # Double sum of probability-weighted absolute differences
        diff_matrix = np.abs(values[:, None] - values[None, :])
        weighted = diff_matrix * (probs[:, None] * probs[None, :])
        numerator = weighted.sum()

        mu = (values * probs).sum()
        return numerator / (2 * mu)

    def compare_normal(self, dist_a, dist_b) -> float:
        mean_a = self.pmf_mean(dist_a)
        mean_b = self.pmf_mean(dist_b)
        std_a = self.pmf_std(dist_a)
        std_b = self.pmf_std(dist_b)

        mean_distance_abs = np.abs(mean_a - mean_b) # type: ignore
        
        variance_sum = std_a + std_b
        
        score = mean_distance_abs/variance_sum
        score = np.clip(score, 0, 1) 

        return score

    def compare_uniform(self, dist_a, dist_b) -> float:
        min_a = self.pmf_min(dist_a)
        min_b = self.pmf_min(dist_b) 
        max_a = self.pmf_max(dist_a)
        max_b = self.pmf_max(dist_b)

        if ((min_a > min_b and min_a > max_b) or (min_b > min_a and min_b > max_a)):
                return 1
        elif ((min_a <= min_b and max_b <= max_a) or (min_b <= min_a and max_a <= max_b)):
                return 0
        elif (min_a > min_b):
            overlap = np.abs(max_a - min_b)
        elif (min_b > min_a):
            overlap = np.abs(max_b - min_a) 
        else:
            overlap = 0.0

        score = (overlap)/((max_a - min_a)+(max_b - min_b))
        score = np.clip(score, 0, 1)
        
        return score
        
    def compare_power(self, dist_a, dist_b) -> float:
        score = np.abs(self.pmf_gini(dist_a) - self.pmf_gini(dist_b))
        score = np.clip(score, 0, 1)
        return score

    def compare_gamma(self, dist_a, dist_b)  -> float:
        alpha1 = (self.pmf_mean(dist_a)**2)/(self.pmf_std(dist_a)**2)
        alpha2 = (self.pmf_mean(dist_b)**2)/(self.pmf_std(dist_b)**2)
        beta1 = (self.pmf_std(dist_a)**2)/(self.pmf_std(dist_a)) #not used
        beta2  = (self.pmf_std(dist_a)**2)/(self.pmf_std(dist_a)) #not used
        score = np.abs(alpha1-alpha2)/10 #10 is hardcoded max difference, higher ones will be treated the same as 10
        score = np.clip(score, 0, 1)
        return score
    
    def stat_dist(self, a:GraphStats, b:GraphStats) -> float:
        a_distribution: dict = getattr(a, self.distributionName)
        b_distribution: dict = getattr(b, self.distributionName)

        dist_type_a = self.identifier.identify_distribution(a_distribution) 
        dist_type_b = self.identifier.identify_distribution(b_distribution) 
        #if distrubutions have diffrent types give max score
        if (dist_type_a != dist_type_b):
            return 1 # -> 1 is the max value (0 is the lowest) 
        elif (dist_type_a == 'normal'):
            return self.compare_normal(a_distribution, b_distribution)
            
        elif (dist_type_a == 'uniform'):
            return self.compare_uniform(a_distribution, b_distribution)   

        elif (dist_type_a == 'powerlaw'):
            return self.compare_power(a_distribution, b_distribution)  

        elif (dist_type_a == 'gamma'): #genaral distribution type - let it be the last resort
            return self.compare_gamma(a_distribution, b_distribution)  
        else :
            return 1.0
