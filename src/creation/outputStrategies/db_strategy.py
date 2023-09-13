import psycopg2
from src.creation.globalSettings import *

def db_strategy(self,algo, user = globalSettings.defaultDBUser, password= globalSettings.defaultDBPassword , host = globalSettings.defaultDBHost, database = globalSettings.defaultDB, **kwargs ):
    df = algo(self, **kwargs) 
    conn = psycopg2.connect(user=user, password=password, host=host, database=database, port = "5432")
    cur = conn.cursor()
    cur.execute(f'''CREATE TABLE IF NOT EXISTS {df.name} 
      (ID BIGSERIAL PRIMARY KEY     NOT NULL,
      r1        INT    NOT NULL,
      r2        INT    NOT NULL,
      uuid      CHAR(36) NOT NULL);''')
    conn.commit()
    cur.execute(f"SELECT uuid FROM {df.name}");
    for idx in df.index:
        cur.execute(f'INSERT INTO {df.name}(r1,r2,uuid) VALUES ({df["r1"][idx]},{df["r2"][idx]},gen_random_uuid())')
        conn.commit()

    conn.close()
    return None



