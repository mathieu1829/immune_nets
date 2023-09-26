import psycopg2
from dotenv import load_dotenv
import os

def db_strategy(self,algo,**kwargs ):
    load_dotenv(); 
    df = algo(self, **kwargs) 
    conn = psycopg2.connect(
        user=os.environ["POSTGRES_USER"],
        password=os.environ["POSTGRES_PASSWORD"],
        host=os.environ["POSTGRES_HOST"],
        database=os.environ["POSTGRES_DB"],
        port = os.environ["POSTGRES_PORT"])
    cur = conn.cursor()
    cur.execute(f'''CREATE TABLE IF NOT EXISTS runs 
      (ID uuid PRIMARY KEY DEFAULT gen_random_uuid(),
      name varchar(255) NOT NULL,
      matrix varchar(255));''')
    conn.commit()
    cur.execute(f'''CREATE TABLE IF NOT EXISTS {df.name} 
      (ID uuid PRIMARY KEY DEFAULT gen_random_uuid(),
      r1        INT    NOT NULL,
      r2        INT    NOT NULL,
      run_id    uuid   NOT NULL,
      CONSTRAINT fk_run_id FOREIGN KEY(run_id) 
      REFERENCES runs(id) ON DELETE CASCADE ON UPDATE CASCADE);''')
    conn.commit()
    matrix = ""
    if "matrix" in kwargs:
      matrix = kwargs["matrix"]
    cur.execute(f'''INSERT INTO runs (name, matrix) VALUES (\'{self.__class__.__name__}\', \'{matrix}\') RETURNING id''')
    uuid = cur.fetchone()[0]
    for idx in df.index:
        cur.execute(f'''INSERT INTO {df.name}(r1,r2,run_id) VALUES ({df["r1"][idx]},{df["r2"][idx]},\'{uuid}')''')
        conn.commit()

    conn.close()
    return None



