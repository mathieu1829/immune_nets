import unittest
from src.orm.db import engine,SessionLocal
from src.orm.models import *
from sqlalchemy import insert,select,delete
from datetime import date,time,datetime
from sqlalchemy.orm import Session

class testORM(unittest.TestCase):
    def test_ORMConnetion(self):
        with engine.connect() as connection:
            pass

    def test_ORMTransaction(self):
        with Session(engine) as session:
            first_record = Dataset(name="a", description="aaa", version="1.1", username="johndoe", modificationDate=datetime.now(), creationDate=datetime.now())
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

