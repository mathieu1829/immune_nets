import pandas as pd

def csv_strategy(self, algo, **kwargs):
    immune_net = algo(self, **kwargs)
    immune_net.to_csv(immune_net.df.name+"_net.csv")
