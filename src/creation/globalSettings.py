import pandas as pd 
from pathlib import Path
from src.creation.io_strategies.df_strategy import df_strategy

class globalSettings(object):
    def __new__(cls):
        if not hasattr(cls, 'instance'):
            cls.instance = super(globalSettings, cls).__new__(cls)
        return cls.instance
    defaultCsvPath = "./"
    defaultDBProvider = "postgres"
    defaultDBUser = "postgres"
    defaultDBPassword = ""
    defaultDBHost = "localhost"
    defaultDB = ""
    defaultOutputStrategy = df_strategy.output
