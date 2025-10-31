import unittest
from src.models import *
from src.db import engine,SessionLocal
from sqlalchemy import insert,select,delete
from datetime import date,time,datetime
from sqlalchemy.orm import Session
from pathlib import Path
import pandas as pd
import uuid

from src.creation.immuneNetwork import ImmuneNetwork
from src.factories import RepertoireFactory, ImmuneNetworkFactory
from src.mappers import RepertoireMapper, NetworkMapper, ImmuneNetworkMapper

class testORM(unittest.TestCase):
    def test_ORMConnetion(self):
        with engine.connect() as connection:
            pass

    def test_ORMSession(self):
        with Session(engine) as session:
            pass

    def test_ORMTransaction(self):
        with Session(engine) as session:
            dataset: Dataset = Dataset(name="test dataset", description="test description")
            session.add(dataset)
            session.commit()

            stmt = select(Dataset).where(Dataset.name == "test dataset")
            result = session.execute(stmt)
            datasetDB: Dataset|None = result.scalars().first()
            self.assertIsNotNone(datasetDB)
            if datasetDB is not None:
                self.assertEqual(datasetDB.description, "test description")
            
            result = session.execute(stmt)
            for obj in result.scalars().all():
                session.delete(obj)
            session.commit()

            result = session.execute(stmt)
            self.assertIsNone(result.scalars().first())
            
    def test_singleRepertoireLoading(self):
        test_dir = Path(__file__).parent / "test_data"
        covid_path = test_dir / "covid_test_clonotypes_0.csv" # covid
        with Session(engine) as session:
            repertoire: Repertoire = RepertoireFactory.fromCSV(name="test_covid",desc="some bile sample", path = covid_path)
            session.add(repertoire)
            session.commit()

            stmt = select(Repertoire).where(Repertoire.name == "test_covid")
            result = session.execute(stmt)
            repertoireDB: Repertoire|None = result.scalars().first()
            self.assertIsNotNone(repertoireDB)
            if repertoireDB is not None:
                self.assertEqual(repertoireDB.description, "some bile sample")

            immuneRepertoire = RepertoireMapper.toImmuneRepertoire(repertoire)
            self.assertEqual(immuneRepertoire.clones.empty,False) 

            result = session.execute(stmt)
            for obj in result.scalars().all():
                session.delete(obj)
            session.commit()

            # result = session.execute(stmt)
            # self.assertIsNone(result.scalars().first())

    # def test_singleNetworkLoading(self):
    #     test_dir = Path(__file__).parent / "test_data"
    #     covid_network_path = test_dir / "covid_test_network_0.pkl" # covid
    #     covid_clonotype_path = test_dir / "covid_test_clonotypes_0.csv" # covid
    #
    #     immuneNetwork: ImmuneNetwork = ImmuneNetworkFactory.fromPickle(covid_network_path)
    #     repertoire: Repertoire = RepertoireFactory.fromCSV(name="test_covid",desc="some bile sample", path = covid_clonotype_path)
    #     repertoire.repertoire_id = uuid.uuid4()
    #     immuneNetwork.sampleId = repertoire.repertoire_id
    #     network: Network = NetworkMapper.fromImmuneNetwork(immuneNetwork)
    #
    #     with Session(engine) as session:
    #         print(f"network id: {network.network_id}")
    #         print(f"network repertoire: {network.repertoire_id}")
    #         session.add(repertoire)
    #         session.add(network)
    #         session.commit()
    #
    #         stmt = select(Network).where(Network.name == "covid test network_0")
    #         stmt2 = select(Repertoire).where(Repertoire.name == "test_covid")
    #         result = session.execute(stmt)
    #         networkDB: Network|None = result.scalars().first()
    #         self.assertIsNotNone(networkDB)
    #         if networkDB is not None:
    #             immuneNetDB: ImmuneNetwork = NetworkMapper.toImmuneNetwork(networkDB)
    #             self.assertEqual(immuneNetwork.graph.empty, False)
    #
    #         result = session.execute(stmt)
    #         for obj in result.scalars().all():
    #             session.delete(obj)
    #         session.commit()
    #         
    #         result = session.execute(stmt2)
    #         for obj in result.scalars().all():
    #             session.delete(obj)
    #         session.commit()
    #
    #         result = session.execute(stmt)
    #         self.assertIsNone(result.scalars().first())
            
if __name__ == '__main__':
    unittest.main()

