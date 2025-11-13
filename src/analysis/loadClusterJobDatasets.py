import argparse
import os

from src.models import Repertoire,Dataset, Network
from src.db import engine
from sqlalchemy.orm import Session, selectinload
from sqlalchemy import select
from pathlib import Path
from src.factories import RepertoireFactory
from src.creation.immuneRepertoire import ImmuneRepertoire
from src.mappers import RepertoireMapper, NetworkMapper
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.distance.alignment import sequenceAligner


def loadClusterJobDatasets(groupPaths):

    immuneRepertoires = {}

    with Session(engine) as session: 
        datasets = {group:Dataset(name=f"{group} dataset", description=" ") for group in groupPaths}
        datasetList = [] 

        for group in groupPaths:
            repertoireList = []
            for file in os.listdir(groupPaths[group]):
                path = groupPaths[group] / file
                metadaGroups = ["group", "id", "description"]
                metadata = { group:data for group, data in zip(metadaGroups, file.split("_"))}
                repertoireList.append(RepertoireFactory.fromCSV(name=f"{metadata['group']} {metadata['id']}",desc=f"{metadata['description']}",path=groupPaths[group]))


            datasets[group].repertoires = repertoireList
            datasetList.append(datasets[group])

        session.add_all(datasetList)
        session.commit()

    print("Loading repertoires: finished")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--group-paths", type=str, required=True)
    args = parser.parse_args()

    groupPaths = args.group_paths
    groupPaths = groupPaths.split(",")
    groupPaths = { pair.split(":")[0]:pair.split(":")[1] for pair in groupPaths}

    loadClusterJobDatasets(groupPaths)

    
    
