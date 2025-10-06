from typing import List
from typing import Optional
from sqlalchemy import ForeignKey
from sqlalchemy import String
from sqlalchemy.orm import DeclarativeBase
from sqlalchemy.orm import Mapped
from sqlalchemy.orm import mapped_column
from sqlalchemy.orm import relationship
import uuid
from sqlalchemy.types import Date
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


# class ClonotypeData(Base): 
#     clone_record_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), 
#                                                   primary_key=True,
#                                                   default=uuid.uuid4)
#     uuid clone_record_id "PK"
#     uuid dataset_id "FK"
#     float proportion "*"
#     string cdr3s_aa "*"
#     string cdr3s_nt "*"
#     string inkt_evidence "*"
#     string mait_evidence "*"
#
# DatasetMetadata {
# uuid dataset_id "PK"
# string name "*"
# string description "*"
# string value "*"
# }
# RepertoireComposition {
# uuid id "PK"
# uuid repertoire_id "FK"
# uuid clone_record_id "FK"
# }
# RepertoireDatasets{
# uuid id "PK"
# uuid repertoire_id "FK"
# uuid dataset_id "FK"
# }
# Repertoire {
# uuid repertoire_id "PK"
# string name "*"
# string description "*"
# string version "*"
# string username "*"
# Date modificationDate "*"
# Date creationDate "*"
# }
# RepertoireStat {
# uuid repertoire_stat_id "PK"
# uuid repertoire_id "FK"
# string stat_name "*"
# string value "*"
# string description "*"
# }
# Network {
# uuid network_id "PK"
# uuid repertoire_id "FK"
# string algorithm "*"
# string distance_function "*"
# string network_algorithm_parameters "*"
# string version "*"
# string username "*"
# Date modificationDate "*"
# Date creationDate "*"
# }
# NetworkStat {
# uuid network_stat_id "PK"
# uuid network_id "FK"
# string stat_name "*"
# string value "*"
# string description "*"
# }
# NetworkData {
# uuid network_record_id "PK"
# uuid network_id "FK"
# int r1 "*"
# int r2 "*"
# string weights "*"
# string r1_weights "*"
# string r2_weights "*"
# }
# Datasets ||--|{ ClonotypeData : "has"
# Datasets ||--|{ DatasetMetadata : "has"
#
# RepertoireComposition }|--|{ ClonotypeData : "has"
# RepertoireDatasets }|--|{ Datasets : "has"
# Repertoire ||--|{ RepertoireComposition : "has"
# Repertoire ||--|{ RepertoireDatasets : "has"
# Repertoire ||--|{ RepertoireStat : "hasRepertoire ||--|{ Network : "creates"
# Network ||--|{ NetworkStat : "has"
# Network ||--|{ NetworkData : "has"

Base.metadata.create_all(engine)
