from sqlalchemy import select
from sqlalchemy.orm import Session, selectinload

from immune_nets.db import engine
from immune_nets.models import Dataset


def get_all_datasets():

    with Session(engine) as session:

        stmt = (
            select(Dataset)
            .options(
                selectinload(Dataset.repertoires)
            )
        )


        return list(
            session.scalars(stmt)
        )

def delete_dataset(dataset_id):
    with Session(engine) as session:

        dataset = session.get(
            Dataset,
            dataset_id
        )

        if not dataset:
            return False

        session.delete(dataset)
        session.commit()
    return True


def create_dataset(name: str,
                   description: str):

    name = name.strip()

    if not name:
        raise ValueError(
            "Dataset name cannot be empty."
        )

    with Session(engine) as session:

        existing = session.scalar(
            select(Dataset)
            .where(Dataset.name == name)
        )

        if existing:
            raise ValueError(
                f"Dataset '{name}' already exists."
            )

        dataset = Dataset(
            name=name,
            description=description
        )

        session.add(dataset)
        session.commit()

def get_dataset(dataset_id):

    with Session(engine) as session:

        stmt = (
            select(Dataset)
            .options(
                selectinload(
                    Dataset.repertoires
                )
            )
            .where(
                Dataset.dataset_id
                == dataset_id
            )
        )

        return session.scalar(stmt)

def get_datasets_by_name(dataset_name):

    with Session(engine) as session:

        stmt = (
            select(Dataset)
            .options(
                selectinload(
                    Dataset.repertoires
                )
            )
            .where(
                Dataset.name.ilike(f"%{dataset_name}%")
            )
        )

        return session.scalar(stmt)
