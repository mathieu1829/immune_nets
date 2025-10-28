import numpy as np
from scipy.stats import wasserstein_distance

from src.analysis.methods.graphStats import GraphStats

class WassersteinStatDistance:
    def __init__(self, distributionName1: str, distributionName2: str):
        self.distributionName1 = distributionName1
        self.distributionName2 = distributionName2
    
    def stat_dist(self, a:GraphStats, b:GraphStats) -> float:
        a_distribution: dict = getattr(a, self.distributionName1)
        b_distribution: dict = getattr(b, self.distributionName2)
        return wasserstein_distance(u_values=list(a_distribution.keys()),
                                    u_weights=list(a_distribution.values()),
                                    v_values=list(b_distribution.keys()),
                                    v_weights=list(b_distribution.values()),
                                   )
        


