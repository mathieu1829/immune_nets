import unittest
from src.orm.models import *
from src.orm.db import engine,SessionLocal
from sqlalchemy import insert,select,delete
from datetime import date,time,datetime
from sqlalchemy.orm import Session
from pathlib import Path
import pandas as pd

class testORM(unittest.TestCase):
    def test_ORMConnetion(self):
        with engine.connect() as connection:
            pass

    def test_ORMTransaction(self):
        with Session(engine) as session:
            first_record = Dataset(name="a", description="aaa")
            session.add(first_record)
            session.commit()

        with engine.connect() as connection:

            stmt = select(Dataset).where(Dataset.name == "a")
            print("After adding:")
            print(list(connection.execute(stmt)))
            stmt = delete(Dataset).where(Dataset.name == "a")
            connection.execute(stmt)
            connection.commit()
            stmt = select(Dataset).where(Dataset.name == "a")
            print("After deleting:")
            print(list(connection.execute(stmt)))
    def test_repertoireLoading(self):
        test_dir = Path(__file__).parent
        covid_path = test_dir / "test_data/covid_test_clonotypes.csv" # covid
        with Session(engine) as session:
            repertoire = Repertoire.from_csv(name="covid",desc="some bile sample", path = covid_path)
            session.add(repertoire)
            session.commit()

        with Session(engine) as session:
            stmt = select(Repertoire).where(Repertoire.name == "covid")
            result = session.execute(stmt)
            repertoire = result.scalars().first()
            print(repertoire.clones)
            results = session.scalars(select(Repertoire).where(Repertoire.name == "covid"))
            for obj in results:
                session.delete(obj)
            session.commit()


            
            

if __name__ == '__main__':
    unittest.main()

