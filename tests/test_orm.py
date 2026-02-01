import unittest
from src.models import *
from src.db import engine,SessionLocal
from sqlalchemy import select
from sqlalchemy.orm import Session
from pathlib import Path

from src.entities import GraphStats
from src.entities import ImmuneNetwork
from src.factories import RepertoireFactory, ImmuneNetworkFactory, ImmuneRepertoireFactory
from src.mappers import RepertoireMapper, NetworkMapper, ImmuneNetworkMapper, NetworkStatMapper

class testORM(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        test_dir = Path(__file__).parent / "test_data"
        cls.covid_network_path = test_dir / "covid_test_network_0.pkl" # covid
        cls.covid_clonotype_path = test_dir / "covid_test_clonotypes_0.csv" # covid
        cls.repertoireName = "test_covid"
        cls.networkName = "covid test network_0"

    def deleteObjects(self, stmt, session):
        result = session.execute(stmt)
        for obj in result.scalars().all():
            session.delete(obj)
        session.commit()

    def confirmDeletion(self, stmt, session):
        result = session.execute(stmt)
        self.assertIsNone(result.scalars().first())

    def deleteAndCheck(self, stmt, session):
        self.deleteObjects(stmt, session)
        self.confirmDeletion(stmt, session)

    def getRepertoire(self) -> Repertoire:
        return RepertoireFactory.fromCSV(name=self.repertoireName,desc="some bile sample", path = self.covid_clonotype_path)

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

            self.deleteAndCheck(stmt, session)
            
    def test_singleRepertoireLoading(self):
        with Session(engine) as session:
            repertoire: Repertoire = self.getRepertoire()
            session.add(repertoire)
            session.commit()

            stmt = select(Repertoire).where(Repertoire.name == self.repertoireName)
            result = session.execute(stmt)
            repertoireDB: Repertoire|None = result.scalars().first()
            self.assertIsNotNone(repertoireDB)
            if repertoireDB is not None:
                self.assertEqual(repertoireDB.description, "some bile sample")

            immuneRepertoire = RepertoireMapper.toImmuneRepertoire(repertoire)
            self.assertEqual(immuneRepertoire.clones.empty,False) 

            self.deleteAndCheck(stmt, session)

    def test_singleNetworkLoading(self):
        immuneNetwork: ImmuneNetwork = ImmuneNetworkFactory.fromPickle(self.covid_network_path)
        
        with Session(engine) as session:
            repertoire = self.getRepertoire()
            session.add(repertoire)
            session.commit()

            immuneNetwork.sampleId = repertoire.repertoire_id
            network: Network = NetworkMapper.fromImmuneNetwork(immuneNetwork)
            session.add(network)
            session.commit()

            stmt = select(Network).where(Network.name == self.networkName)
            stmt2 = select(Repertoire).where(Repertoire.name == self.repertoireName)

            result = session.execute(stmt)
            networkDB: Network|None = result.scalars().first()
            self.assertIsNotNone(networkDB)
            if networkDB is not None:
                immuneNetDB: ImmuneNetwork = NetworkMapper.toImmuneNetwork(networkDB)
                self.assertEqual(immuneNetDB.graph.empty, False)

            self.deleteAndCheck(stmt, session)
            self.deleteAndCheck(stmt2, session)

    def test_cascadeRepertoireDelete(self):

        
        immuneNetwork: ImmuneNetwork = ImmuneNetworkFactory.fromPickle(self.covid_network_path)
        
        with Session(engine) as session:
            repertoire = self.getRepertoire()

            immuneNetwork.sampleId = repertoire.repertoire_id
            network: Network = NetworkMapper.fromImmuneNetwork(immuneNetwork)
            repertoire.repertoire_networks.append(network)
            session.add(repertoire)
            session.commit()

            stmt = select(Network).where(Network.name == self.networkName)
            stmt2 = select(Repertoire).where(Repertoire.name == self.repertoireName)

            result = session.execute(stmt)
            networkDB: Network|None = result.scalars().first()
            self.assertIsNotNone(networkDB)
            if networkDB is not None:
                immuneNetDB: ImmuneNetwork = NetworkMapper.toImmuneNetwork(networkDB)
                self.assertEqual(immuneNetDB.graph.empty, False)

            self.deleteObjects(stmt2, session)

            self.confirmDeletion(stmt, session)
            self.confirmDeletion(stmt2, session)

    def test_networkStatLoading(self):
        immuneNetwork: ImmuneNetwork = ImmuneNetworkFactory.fromPickle(self.covid_network_path)
        stats = GraphStats(immuneNetwork)
        
        with Session(engine) as session:
            repertoire = self.getRepertoire()

            immuneNetwork.sampleId = repertoire.repertoire_id
            network: Network = NetworkMapper.fromImmuneNetwork(immuneNetwork)
            network.network_stats = NetworkStatMapper.fromGraphStat(stats,immuneNetwork)
            repertoire.repertoire_networks.append(network)
            session.add(repertoire)
            session.commit()

            stmt = select(Network).where(Network.name == self.networkName)
            stmt2 = select(Repertoire).where(Repertoire.name == self.repertoireName)
            stmt3 = select(NetworkStat).where(NetworkStat.network_id == network.network_id)

            result = session.execute(stmt)
            networkDB: Network|None = result.scalars().first()
            self.assertIsNotNone(networkDB)
            if networkDB is not None:
                statsDB = NetworkStatMapper.toGraphStat(networkDB.network_stats)
                self.assertEqual(statsDB.toList(), stats.toList())

            self.deleteObjects(stmt2, session)

            self.confirmDeletion(stmt, session)
            self.confirmDeletion(stmt2, session)
            self.confirmDeletion(stmt3, session)
            
if __name__ == '__main__':
    unittest.main()

