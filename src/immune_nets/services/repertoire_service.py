from sqlalchemy import select
from sqlalchemy.orm import Session, selectinload

from immune_nets.db import engine
from immune_nets.models import Dataset, Repertoire
from immune_nets.factories import RepertoireFactory
import tempfile

def add_repertoire_to_dataset(
        dataset_id,
        name,
        desc,
        uploaded_file):
    with Session(engine) as session:

        dataset = session.get(
            Dataset,
            dataset_id
        )

        if dataset is None:
            raise ValueError(
                "Dataset does not exist."
            )

        with tempfile.NamedTemporaryFile(
                suffix=".csv",
                delete=False
        ) as tmp:

            uploaded_file.save(
                tmp.name
            )

            repertoire = (
                RepertoireFactory
                .fromCSV(
                    name=name,
                    desc=desc,
                    path=tmp.name
                )
            )

        dataset.repertoires.append(
            repertoire
        )

        session.add(repertoire)

        session.commit()

def delete_repertoire(repertoire_id):
    with Session(engine) as session:

        repertoire = session.get(
            Repertoire,
            repertoire_id
        )

        if not repertoire:
            return False

        session.delete(repertoire)
        session.commit()
    return True

def get_repertoire(repertoire_id):

    with Session(engine) as session:

        stmt = (
            select(Repertoire)
            .options(
                selectinload(
                    Repertoire.clonotypes
                )
            )
            .where(
                Repertoire.repertoire_id
                == repertoire_id
            )
        )

        return session.scalar(stmt)

def get_all_repertoires():

    with Session(engine) as session:

        stmt = (
            select(Repertoire)
            .options(
                selectinload(
                    Repertoire.clonotypes
                )
            )
        )


        return list(
            session.scalars(stmt)
        )

