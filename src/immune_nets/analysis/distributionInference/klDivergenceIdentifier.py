import numpy as np

from .distributionIdentifier import DistributionIdentifier

class KlDivergenceIdentifier(DistributionIdentifier):
    """
    Identifies distributions using Kullback–Leibler divergence.

    Generates samples from the input distribution and candidate theoretical
    distributions, then compares them using KL-divergence.
    """

    def kl_div(self, P, Q, eps=1e-12):
        """
        Calculate the KL-divergence between two discrete distributions.

        :param P: Empirical probability distribution.
        :param Q: Theoretical probability distribution.
        :param eps: Small value to avoid division by zero.
        :return: KL-divergence value.
        """
        P = P + eps
        Q = Q + eps
        P /= P.sum()
        Q /= Q.sum()
        return np.sum(P * np.log(P / Q))

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

                results[name] = self.kl_div(probs, theoretical)
                # if that fails assign high score to guarantee choosing different distribution
            except:
                results[name] = 10
        best = min(results, key=results.get) # type: ignore
        
        return best
