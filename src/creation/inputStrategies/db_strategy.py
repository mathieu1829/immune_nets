from dotenv import load_dotenv
import os
import psycopg2
import pandas as pd
from src.creation.algorithms.common_methods import *

def db_strategy(name):
    load_dotenv(); 
    conn = psycopg2.connect(
        user=os.environ["POSTGRES_USER"],
        password=os.environ["POSTGRES_PASSWORD"],
        host=os.environ["POSTGRES_HOST"],
        database=os.environ["POSTGRES_DB"],
        port = os.environ["POSTGRES_PORT"])
    cur  =  conn.cursor()
    cur.execute(f'''SELECT clonotype_id, frequency, proportion, cdr3s_aa, cdr3s_nt, inkt_evidence, mait_evidence FROM {name}''')
    data = cur.fetchall()
    row_list = []
    for row in data:
        row_list.append(row)
    df = pd.DataFrame(row_list, columns=["clonotype_id", "frequency", "proportion", "cdr3s_aa", "cdr3s_nt", "inkt_evidence", "mait_evidence"])
    cur.close()
    df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
    df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))
    df.name = "networks"
    return df