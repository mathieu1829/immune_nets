from scipy.spatial.distance import hamming


class negativeLevenshteinDistance:
    def tcr_dist(self,x,ref):
        maxLen = max(len(x), len(ref))
        x = x.ljust(maxLen)
        ref = ref.ljust(maxLen)
        return -sum([ i != j for i,j in zip(x,ref)])

    def __str__(self):
        return "negative_levenshtein"
