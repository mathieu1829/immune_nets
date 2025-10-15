
from src.models import Repertoire,Dataset, Network
from src.db import engine
from sqlalchemy.orm import Session, selectinload
from sqlalchemy import select
from pathlib import Path
from src.factories import ImmuneRepertoireFactory
from src.creation.immuneRepertoire import ImmuneRepertoire
from src.mappers import RepertoireMapper, NetworkMapper
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.distance.alignment import sequenceAligner


root_dir = Path(__file__).parent.parent.parent

leukemia_path = root_dir  / "tests/test_data/leukemia_test_clonotypes.csv" # leukemia
covid_path = root_dir / "tests/test_data/covid_test_clonotypes.csv" # covid
healthy_path = root_dir / "tests/test_data/healthy_test_clonotypes_1.csv" #healthy

groups = ["healthy", "leukemia", "covid"]
paths = [healthy_path, leukemia_path, covid_path]
groupPaths = {group:path for group,path in zip(groups,paths)}

immuneRepertoires = {}

with Session(engine) as session: 
    datasets = {group:Dataset(name=f"{group} test dataset", description=" ") for group in groups}
    datasetList = [] 

    for group in groups:
        name = f"{group} test dataset"
        stmt = select(Dataset).options(selectinload(Dataset.repertoires).selectinload(Repertoire.clonotypes)).where(Dataset.name == name)
        result = session.execute(stmt)
        dataset = result.scalars().first()

        if dataset is None:
            print(f"\tDataset: {name} - not present the db; initializing ...")
            immuneRepertoires[group] = ImmuneRepertoireFactory.fromCSV(name=f"{group} test repertoire",desc=" ",path=groupPaths[group])
            datasets[group].repertoires = [
                RepertoireMapper.fromImmuneRepertoire(immuneRepertoires[group])
                ]
            datasetList.append(datasets[group])
        else: 
            print(f"\tDataset: {name} - is present the db; skipping generation ...")
            datasetList.append(dataset)
            immuneRepertoires[group] = RepertoireMapper.toImmuneRepertoire(dataset.repertoires[0])

    session.add_all(datasetList)
    session.commit()

print("Loading repertoires: finished")
print("Computing test networks:")

distance_fun = sequenceAligner("BLOSUM62")
immuneNets = {}

with Session(engine) as session:
    networks = []
    for group in groups:
        name = f"{group} test network"
        stmt = select(Network).where(Network.name == name)
        result = session.execute(stmt)
        network = result.scalars().first()
        if network is None:
            print(f"\tNetwork: {name} - not present the db; initializing ...")
            net = simpleBetaDistance(repertoire=immuneRepertoires[group],
                                        distance=distance_fun,
                                        threshold=0.3)
            net.name = name
            immuneNets[group] = net
            networks.append(NetworkMapper.fromImmuneNetwork(immuneNets[group]))

            print(f"\tNetwork: {name} - computed")
        else :
            print(f"\tNetwork: {name} - is present the db; skipping generation ...")
            networks.append(network)

    session.add_all(networks)
    session.commit()
    print("Loading networks: finished")
