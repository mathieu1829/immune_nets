from dotenv import load_dotenv
import os
import psycopg2
import pandas as pd
from src.creation.algorithms.common_methods import *


load_dotenv() 
class db_strategy():
    username = os.environ["POSTGRES_USER"]
    password = os.environ["POSTGRES_PASSWORD"]
    database = os.environ["POSTGRES_DB"]
    host = os.environ["POSTGRES_HOST"]
    port = os.environ["POSTGRES_PORT"]
    def input(self):
        conn = psycopg2.connect(user=self.username,password=self.password, port=self.password, host=self.host, database=self.database)
        cur  =  conn.cursor()
        cur.execute(f'''SELECT clonotype_id, frequency, proportion, cdr3s_aa, cdr3s_nt, inkt_evidence, mait_evidence FROM data''')
        data = cur.fetchall()
        row_list = []
        for row in data:
            row_list.append(row)
        df = pd.DataFrame(row_list, columns=["clonotype_id", "frequency", "proportion", "cdr3s_aa", "cdr3s_nt", "inkt_evidence", "mait_evidence"])
        cur.close()
        df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
        df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
        df.name = "data"
        return df
    @staticmethod
    def output(self, algo, **kwargs):
        load_dotenv(); 
        username = os.environ["POSTGRES_USER"]
        password = os.environ["POSTGRES_PASSWORD"]
        database = os.environ["POSTGRES_DB"]
        host = os.environ["POSTGRES_HOST"]
        df = algo(self, **kwargs) 
        conn = psycopg2.connect(user=username,password=password, port=password, host=host, database=database)
        cur = conn.cursor()
        cur.execute(f'''CREATE TABLE IF NOT EXISTS runs 
        (ID uuid PRIMARY KEY DEFAULT gen_random_uuid(),
        name varchar(255) NOT NULL,
        method varchar(255) NOT NULL,
        distance varchar(255));''')
        conn.commit()
        cur.execute(f'''CREATE TABLE IF NOT EXISTS results 
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
        
        cur.execute(f'''INSERT INTO runs (name, method, distance) VALUES (\'{df.name}\', \'{self.__class__.__name__}\', \'{matrix}\') RETURNING id''')
        # In case we decided to name runs with 
        # cur.execute(f'''INSERT INTO runs (name, matrix) VALUES (\'{self.__class__.__name__}\', \'{matrix}\') RETURNING id''')
        uuid = cur.fetchone()[0]
        for idx in df.index:
            cur.execute(f'''INSERT INTO results(r1,r2,run_id) VALUES ({df["r1"][idx]},{df["r2"][idx]},\'{uuid}')''')
            conn.commit()

        conn.close()
        return None
