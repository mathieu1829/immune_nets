import pandas as pd
import numpy as np
import random

class testGenerator:
    def generate_random_cdr3(self,min_length = 11,max_length=16):
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        random_length = random.randint(min_length, max_length)
        return ''.join(random.choices(amino_acids, k=random_length))

    def mutateCDR3(self, cdr3seq, numMutations=1):
        seq = list(cdr3seq)
        for _ in range(numMutations):
            mutationIdx = random.randint(0,len(seq)-1)
            seq[mutationIdx] = random.choice("ACDEFGHIKLMNPQRSTVWY".replace(seq[mutationIdx],""))
        return "".join(seq)


if __name__ == "__main__":
    testGen = testGenerator()
    randomCDR3 = testGen.generate_random_cdr3()
    print(f"random cdr3: {randomCDR3}")
    for _ in range(10):
        print(f"mutated random cdr3: {testGen.mutateCDR3(randomCDR3)}")
    
