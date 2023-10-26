import pandas as pd

class df_strategy:
    def output(self, algo, **kwargs):
        return algo(self,**kwargs) #if it doesn't work the arguments can be specified manually
