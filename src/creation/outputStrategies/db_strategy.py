import pandas as pd
import psycopg2
import uuid
import numpy as np
from src.creation.globalSettings import *

def db_strategy(self,algo, user = globalSettings.defaultDBUser, password= globalSettings.defaultDBPassword , host = globalSettings.defaultDBHost, database = globalSettings.defaultDB, **kwargs ):
    df = algo(self, **kwargs) 
    conn = psycopg2.connect(user=user, password=password, host=host, database=database, port = "5432")
    cur = conn.cursor()
    cur.execute(f'''CREATE TABLE IF NOT EXISTS {df.name} 
      (ID INT PRIMARY KEY     NOT NULL,
      r1        INT    NOT NULL,
      r2        INT    NOT NULL,
      uuid      CHAR(32)    NOT NULL);''')
    conn.commit()
    nextUUID = uuid.uuid4().hex 
    cur.execute(f"SELECT uuid FROM {df.name}");
    UUIDS = np.unique(np.array([x[0] for x in cur.fetchall()]))
    nextID = cur.rowcount 
    while(len(UUIDS) and (UUIDS == nextUUID).any()):
        nextUUID = uuid.uuid4().hex
    for idx in df.index:
        nextID+=1
        cur.execute(f'INSERT INTO {df.name}(ID,r1,r2,uuid) VALUES ({nextID},{df["r1"][idx]},{df["r2"][idx]},\'{nextUUID}\')')
        conn.commit()

    conn.close()
    return None



