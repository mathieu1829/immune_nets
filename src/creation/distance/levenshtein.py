

class levenshteinDistance:
    def tcr_dist(self,x,ref):
        maxLen = max(len(x), len(ref))
        x = x.ljust(maxLen)
        ref = ref.ljust(maxLen)
        return sum([ i != j for i,j in zip(x,ref)])
        # return -hamming(list(x), list(ref))

    def __str__(self):
        return "levenshtein"

if __name__ == "__main__":
    a = levenshteinDistance()
    x = "aa"
    ref = "aaa"
    print(a.tcr_dist(x,ref))
