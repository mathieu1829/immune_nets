import numpy as np
from scipy import stats as scipy_stats

from .distributionIdentifier import DistributionIdentifier

class KolomogorovIndentifier(DistributionIdentifier):
    """
    Identifies distributions using Kolomogorov-Smirov test.

    Generates samples from the input distribution, fits :candidate theoretical
    distributions, then compares them using either Kolomogorov-Smirov value or pvalue.

    """
    def __init__(self, usePvalue=True):
        """
        :param usePvalue: If True, compare distributions using p-value; otherwise use the KS statistic. 
        """
        self.usePvalueFlag = usePvalue

    def identify_distribution(self, statDict):
        values = np.array(list(statDict.keys()))
        probs  = np.array(list(statDict.values()))
        if probs.sum() != 1:
            probs = probs / probs.sum()

        results = {}
        for name, dist in self.candidates.items():
            # generate samples from discrete distribution
            sample = np.random.choice(values, size=100000, p=probs)
            # try to fit the to theoretical candidate distribution
            try:
                params = dist.fit(sample)

                theoretical = dist.pdf(values, *params)
                theoretical /= theoretical.sum()
               
                statistic, pvalue = scipy_stats.kstest(sample, name, params)

                results[name] = pvalue if self.usePvalueFlag else statistic
                # if that fails assign high score to guarantee choosing different distribution
            except:
                results[name] = 100000
        if self.usePvalueFlag :
            best = max(results, key=results.get) # type: ignore
        else:
            best = min(results, key=results.get) # type: ignore

        
        return best
