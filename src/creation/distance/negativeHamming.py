from scipy.spatial.distance import hamming


class negativeHammingDistance:
    def tcr_dist(self,x,ref):
        longer,shorter = (x,ref) if len(x) > len(ref) else (ref,x)
        shorter ="".join([shorter,"".join([" " for i in range(len(longer)-len(shorter))])])
        return sum([ i != j for i,j in zip(longer,shorter)])
        # return hamming(list(longer), list(shorter))
    def __str__(self):
        return "negative_hamming"
