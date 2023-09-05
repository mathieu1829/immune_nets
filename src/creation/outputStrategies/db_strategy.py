import pandas as pd
from pony.orm import *
from src.creation.globalSettings import *

# def db_strategy(self,algo, user = globalSettings.defaultDBUser, password= globalSettings.defaultDBPassword , host = globalSettings.defaultDBHost, database = globalSettings.defaultDB, **kwargs ):

@db_session
def db_strategy(self,algo, **kwargs ):
    df = algo(self, **kwargs) 
    db = Database()
    # db.bind(provider='postgres', user=user, password=password, host=host, database=database)
    # exec(f'class {df.name}(db.Entity):\n\tr1=Required(int)\n\tr2=Required(int)') 
    # db.generate_mapping(create_tables=True)
    for idx in df.index:
        exec(f'db_record = {df.name}(r1 = {df["r1"][idx]}, r2 = {df["r2"][idx]})')
        commit()
    return None



