from typing import List
from typing import Optional
from sqlalchemy import ForeignKey
from sqlalchemy.orm import DeclarativeBase
from sqlalchemy.orm import Mapped
from sqlalchemy.orm import mapped_column
from sqlalchemy.orm import relationship
import uuid
from sqlalchemy.types import Date, Float, String, Integer
from sqlalchemy.dialects.postgresql import UUID
from datetime import date,time
from src.orm.db import engine

class Base(DeclarativeBase):
     pass

class Dataset(Base):
    __tablename__ = "dataset"

    dataset_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                  primary_key=True,
                                                  default=uuid.uuid4)
    name: Mapped[str] =  mapped_column(String(30), nullable = False)
    description: Mapped[str] =  mapped_column(String(120), nullable = False)
    version: Mapped[str] =  mapped_column(String(10), nullable = False)
    username: Mapped[str] =  mapped_column(String(30), nullable = False)
    modificationDate: Mapped[date] =  mapped_column(Date, nullable = False)
    creationDate: Mapped[date] =  mapped_column(Date, nullable = False)

    dataset_metadata_list: Mapped[List["DatasetMetadata"]] = relationship(back_populates="source_dataset")
    repertoires: Mapped[List["Repertoire"]] = relationship(
        secondary="repertoire_datasets",   # table name as string
        back_populates="datasets",
    )




class ClonotypeData(Base): 
    __tablename__ = "clonotype_data"

    clone_record_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                            primary_key=True,
                                                            default=uuid.uuid4)
    dataset_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("dataset.dataset_id"))
    cdr3s_aa: Mapped[str] = mapped_column(String(100))
    cdr3s_nt: Mapped[str] = mapped_column(String(300))
    proportion: Mapped[float] = mapped_column(Float)
    inkt_evidence: Mapped[str] = mapped_column(String(30))
    mait_evidence: Mapped[str] = mapped_column(String(30))

    source_repertoires: Mapped[List["Repertoire"]] = relationship(
        secondary="repertoire_composition",   # table name as string
        back_populates="clonotypes",
    )

class DatasetMetadata(Base): 
    __tablename__ = "dataset_metadata"

    dataset_metadata_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                           primary_key=True,
                                                           default=uuid.uuid4)
    dataset_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("dataset.dataset_id"))
    name: Mapped[str] =  mapped_column(String(30), nullable = False)
    description: Mapped[str] =  mapped_column(String(120), nullable = False)
    value: Mapped[str] =  mapped_column(String(120), nullable = False)


    source_dataset: Mapped["Dataset"] = relationship(back_populates="dataset_metadata_list")


class RepertoireComposition(Base):
    __tablename__ = "repertoire_composition"

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                          primary_key=True,
                                          default=uuid.uuid4)
    repertoire_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("repertoire.repertoire_id"))
    clone_record_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("clonotype_data.clone_record_id"))


class RepertoireDatasets(Base):
    __tablename__ = "repertoire_datasets"

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                          primary_key=True,
                                          default=uuid.uuid4)
    repertoire_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("repertoire.repertoire_id"))
    dataset_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("dataset.dataset_id"))

class Repertoire(Base):
    __tablename__ = "repertoire"

    repertoire_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                     primary_key=True,
                                                     default=uuid.uuid4)
    name: Mapped[str] =  mapped_column(String(30), nullable = False)
    description: Mapped[str] =  mapped_column(String(120), nullable = False)
    version: Mapped[str] =  mapped_column(String(10), nullable = False)
    username: Mapped[str] =  mapped_column(String(30), nullable = False)
    modificationDate: Mapped[date] =  mapped_column(Date, nullable = False)
    creationDate: Mapped[date] =  mapped_column(Date, nullable = False)

    datasets: Mapped[List["Dataset"]] = relationship(
        secondary="repertoire_datasets",   # table name as string
        back_populates="repertoires",
    )
    clonotypes: Mapped[List["ClonotypeData"]] = relationship(
        secondary="repertoire_composition",   # table name as string
        back_populates="source_repertoires",
    )
    repertoire_stats: Mapped[List["RepertoireStat"]] = relationship(back_populates="source_repertoire")
    repertoire_networks: Mapped[List["Network"]] = relationship(back_populates="source_repertoire")

class RepertoireStat(Base):
    __tablename__ = "repertoire_stat_id"

    repertoire_stat_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                          primary_key=True,
                                                          default=uuid.uuid4)
    repertoire_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("repertoire.repertoire_id"))
    stat_name: Mapped[str] =  mapped_column(String(30), nullable = False)
    value: Mapped[str] =  mapped_column(String(120), nullable = False)
    description: Mapped[str] =  mapped_column(String(120), nullable = False)

    source_repertoire: Mapped["Repertoire"] = relationship(back_populates="repertoire_stats")

class Network(Base):
    __tablename__ = "network"

    network_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                          primary_key=True,
                                                          default=uuid.uuid4)
    repertoire_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("repertoire.repertoire_id"))

    algorithm: Mapped[str] =  mapped_column(String(30), nullable = False)
    distance_function: Mapped[str] =  mapped_column(String(30), nullable = False)
    network_algorithm_parameters: Mapped[str] =  mapped_column(String(30), nullable = False)
    version: Mapped[str] =  mapped_column(String(10), nullable = False)
    username: Mapped[str] =  mapped_column(String(30), nullable = False)
    modificationDate: Mapped[date] =  mapped_column(Date, nullable = False)
    creationDate: Mapped[date] =  mapped_column(Date, nullable = False)

    source_repertoire: Mapped["Repertoire"] = relationship(back_populates="repertoire_networks")
    network_stats: Mapped[List["NetworkStat"]] = relationship(back_populates="source_network")
    network_edges: Mapped[List["NetworkData"]] = relationship(back_populates="source_network")

class NetworkStat(Base):
    __tablename__ = "network_stat"

    network_stat_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                          primary_key=True,
                                                          default=uuid.uuid4)
    network_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("network.network_id"))
    stat_name: Mapped[str] =  mapped_column(String(30), nullable = False)
    description: Mapped[str] =  mapped_column(String(120), nullable = False)
    value: Mapped[str] =  mapped_column(String(120), nullable = False)

    source_network: Mapped["Network"] = relationship(back_populates="network_stats")

class NetworkData(Base):
    __tablename__ = "network_data"

    network_record_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
                                                          primary_key=True,
                                                          default=uuid.uuid4)
    network_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("network.network_id"))
    r1: Mapped[int] = mapped_column(Integer)
    r2: Mapped[int] = mapped_column(Integer)
    weights: Mapped[str] = mapped_column(String(50))
    r1_weights: Mapped[str] = mapped_column(String(50))
    r2_weights: Mapped[str] = mapped_column(String(50))


    source_network: Mapped["Network"] = relationship(back_populates="network_edges")


# Datasets ||--|{ ClonotypeData : "has"
# Datasets ||--|{ DatasetMetadata : "has"
#
# RepertoireComposition }|--|{ ClonotypeData : "has"
# RepertoireDatasets }|--|{ Datasets : "has"
# Repertoire ||--|{ RepertoireComposition : "has"
# Repertoire ||--|{ RepertoireDatasets : "has"

# Repertoire ||--|{ RepertoireStat : "hasRepertoire
# Repertoire ||--|{ Network : "creates"
# Network ||--|{ NetworkStat : "has"
# Network ||--|{ NetworkData : "has"

Base.metadata.create_all(engine)
