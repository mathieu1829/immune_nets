from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from src.orm.models import Repertoire,Dataset,Network
from src.orm.db import engine
from sqlalchemy.orm import Session, selectinload
from sqlalchemy import insert,select,delete
from src.creation.immuneRepertoire import immuneRepertoire
from src.analysis.methods.graphletComposition import graphletComposition
from src.creation.distance.alignment import sequenceAligner
import time


with Session(engine) as session: 
    stmt1 = select(Dataset).options(selectinload(Dataset.repertoires).selectinload(Repertoire.clonotypes)).where(Dataset.name == "healthy")
    result1 = session.execute(stmt1)
    # print(list(result1.all()))
    healthy_dataset = result1.scalars().first()
    db_repertoire = healthy_dataset.repertoires[0]
    normal_repertoire = immuneRepertoire(db_repertoire)


if __name__ == '__main__':
    distanceFun = sequenceAligner("BLOSUM62")
    
    start = time.time()
    network = simpleBetaDistance(repertoire=db_repertoire, distance=distanceFun)
    end = time.time()
    print(f"network simple db execution time: {end - start:.4f} seconds")

    start = time.time()
    stats = graphletComposition(network)
    end = time.time()
    print(f"stats simple db execution time: {end - start:.4f} seconds")
    print()

    start = time.time()
    network = simpleVectorBetaDistance(repertoire=db_repertoire, distance=distanceFun)
    end = time.time()
    print(f"network vector db execution time: {end - start:.4f} seconds")

    start = time.time()
    stats = graphletComposition(network)
    end = time.time()
    print(f"stats vector db execution time: {end - start:.4f} seconds")
    print()

    start = time.time()
    network = simpleBetaDistance(repertoire=normal_repertoire, distance=distanceFun, alt_rep=True)
    end = time.time()
    print(f"network simple normal execution time: {end - start:.4f} seconds")

    start = time.time()
    stats = graphletComposition(network)
    end = time.time()
    print(f"stats simple normal execution time: {end - start:.4f} seconds")
    print()

    start = time.time()
    network_normal = simpleVectorBetaDistance(repertoire=normal_repertoire, distance=distanceFun, alt_rep=True)
    end = time.time()
    print(f"network vector normal execution time: {end - start:.4f} seconds")

    start = time.time()
    stats_normal = graphletComposition(network_normal)
    end = time.time()
    print(f"stats vector normal execution time: {end - start:.4f} seconds")
    print()



