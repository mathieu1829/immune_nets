import pandas as pd

class df_strategy:
    def output(self, algo, **kwargs):
        return algo(**kwargs) #if it doesn't work the arguments can be specified manually
