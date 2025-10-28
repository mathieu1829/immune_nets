
from src.models import Repertoire,Dataset, Network
from src.db import engine
from sqlalchemy.orm import Session, selectinload
from sqlalchemy import select
from pathlib import Path
from src.factories import ImmuneRepertoireFactory, ImmuneNetworkFactory
from src.creation.immuneRepertoire import ImmuneRepertoire
from src.mappers import RepertoireMapper, NetworkMapper
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.distance.alignment import sequenceAligner


root_dir = Path(__file__).parent.parent.parent
test_data_path = root_dir / "tests/test_data"

leukemia_path = test_data_path  / "leukemia_test_clonotypes_0.csv" # leukemia
covid_path = test_data_path / "covid_test_clonotypes_0.csv" # covid
healthy_path = test_data_path / "healthy_test_clonotypes_1.csv" #healthy

leukemia_network_path = test_data_path  / "leukemia_test_network.pkl" # leukemia
covid_network_path = test_data_path / "covid_test_network.pkl" # covid
healthy_network_path = test_data_path / "healthy_test_network.pkl" #healthy

groups = ["healthy", "leukemia", "covid"]
paths = [healthy_path, leukemia_path, covid_path]
networkPaths = [healthy_network_path, leukemia_network_path, covid_network_path]
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
    for group, networkPath in zip(groups, networkPaths):
        name = f"{group} test network"
        stmt = select(Network).where(Network.name == name)
        result = session.execute(stmt)
        network = result.scalars().first()
        newNet = None

        if network is not None :
            print(f"\tNetwork: {name} - is present the db; skipping generation ...")
            networks.append(network)
            continue

        if newNet is None and networkPath.is_file():
            print(f"\tNetwork: {name} - not present the db; retrieving from file ...")
            newNet = ImmuneNetworkFactory.fromPickle(networkPath)
            newNet.sampleId = immuneRepertoires[group].repertoire_id

        if newNet is None:
            print(f"\tNetwork: {name} - not present the db; initializing ...")
            newNet = simpleBetaDistance(repertoire=immuneRepertoires[group],
                                        distance=distance_fun,
                                        threshold=0.3)
            newNet.name = name 

        networks.append(NetworkMapper.fromImmuneNetwork(newNet))

        print(f"\tNetwork: {name} - computed")


    session.add_all(networks)
    session.commit()
    print("Loading networks: finished")
