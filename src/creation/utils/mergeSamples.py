import pandas as pd

def mergeSamples(samples):
    res = samples[0]
    for sample in samples:
        if res.equals(sample):
            continue;
        else:
            res.append(sample);
