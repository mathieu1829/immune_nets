from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from pathlib import Path
from dotenv import load_dotenv
import os 
from immune_nets.models import *



load_dotenv()
pgSocket = os.environ["PGHOST"]
pg_db = os.environ["POSTGRES_DB"]
pg_user = os.environ["POSTGRES_USER"]
pg_passwd = os.environ["POSTGRES_PASSWORD"]
pg_passwd = f":{pg_passwd}" if pg_passwd != "" else pg_passwd

rootDir = Path(__file__).parent.parent.parent
pgSocketPath = rootDir / pgSocket


engine = create_engine(f"postgresql+psycopg2://{pg_user}{pg_passwd}@/{pg_db}?host={pgSocketPath}")
SessionLocal = sessionmaker(bind=engine, autoflush=False, autocommit=False)



Base.metadata.create_all(engine)
