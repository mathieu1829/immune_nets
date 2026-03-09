import argparse
import os

from src.models import Repertoire,Dataset, Network
from src.db import engine
from sqlalchemy.orm import Session, selectinload
from sqlalchemy import select
from pathlib import Path
from src.factories import RepertoireFactory, ImmuneRepertoireFactory
from src.entities import ImmuneRepertoire
from src.mappers import RepertoireMapper, NetworkMapper
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.distance.alignment import sequenceAligner
from src.cluster.utils.commonMethods import loadDatasetsFromFile

#move to some utils
def max_string_lengths_across_dfs(repertoires):
    # Find string / object columns (assume same across all dfs)
    string_cols = repertoires[0].clones.select_dtypes(include=['object', 'string']).columns
    
    max_lengths = {col: 0 for col in string_cols}

    for repertoire in repertoires:
        for col in string_cols:
            col_max = repertoire.clones[col].astype(str).apply(len).max()
            if col_max > max_lengths[col]:
                max_lengths[col] = col_max

    return max_lengths


def loadClusterJobDatasetsToDB(groupPaths):

    immuneRepertoires = {}

    with Session(engine) as session: 
        datasets = {group:Dataset(name=f"{group} dataset", description=" ") for group in groupPaths}
        datasetList = [] 
        repertoireDatasets = loadDatasetsFromFile(groupPaths, db=True) 

        for group in repertoireDatasets:
            # print(max_string_lengths_across_dfs(repertoireList))
            datasets[group].repertoires = repertoireDatasets[group]
            datasetList.append(datasets[group])

        session.add_all(datasetList)
        session.commit()

        print("Loading repertoires: finished")
        print()

def showDBContents(groupPaths):
    with Session(engine) as session: 
        print("Database contents:")
        for group in groupPaths:
            print(f"\t{group} dataset:")
            stmt = select(Dataset).options(selectinload(Dataset.repertoires).selectinload(Repertoire.clonotypes)).where(Dataset.name == f"{group} dataset")
            result = session.execute(stmt)
            dataset = result.scalars().first()

            for rep in dataset.repertoires:
                print(f"\t\t{rep.name}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--group-paths", type=str, required=True)
    args = parser.parse_args()

    groupPaths = args.group_paths
    groupPaths = groupPaths.split(",")
    groupPaths = { pair.split(":")[0]:pair.split(":")[1] for pair in groupPaths}

    loadClusterJobDatasetsToDB(groupPaths)
    showDBContents(groupPaths)



    
    
