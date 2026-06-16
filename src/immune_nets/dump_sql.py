from sqlalchemy import create_engine
# Import your Base and all models so they register with the metadata
from immune_nets.models.base import Base
from immune_nets.models.dataset import Dataset
from immune_nets.models.repertoire import Repertoire
from immune_nets.models.clonotypeData import ClonotypeData
from immune_nets.models.network import Network
from immune_nets.models.networkData import NetworkData
from immune_nets.models.networkStat import NetworkStat
from immune_nets.models.repertoireDatasets import RepertoireDatasets
from immune_nets.models.repertoireStat import RepertoireStat

def dump_ddl():
    # Use a dummy postgresql URL but with a custom strategy to catch the SQL
    engine = create_engine("postgresql://", strategy="mock", executor=lambda sql, *multiparams, **params: print(f"{sql};"))
    
    print("-- Automated DB Creation Script generated from ORM Models\n")
    Base.metadata.create_all(engine)

if __name__ == "__main__":
    dump_ddl()
