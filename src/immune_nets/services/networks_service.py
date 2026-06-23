from sqlalchemy import select
from sqlalchemy.orm import Session, selectinload

from immune_nets.db import engine
from immune_nets.models import Network, Repertoire

from immune_nets.mappers import NetworkMapper

def get_all_networks():

    with Session(engine) as session:

        stmt = (
            select(Network)
            .options(
                selectinload(Network.network_edges),
                selectinload(Network.source_repertoire)
            )
        )


        return list(
            session.scalars(stmt)
        )

def delete_network(network_id):
    with Session(engine) as session:

        network = session.get(
            Network,
            network_id 
        )

        if not network:
            return False

        session.delete(network)
        session.commit()
    return True

def add_network(immuneNetwork):
    with Session(engine) as session:
        network = NetworkMapper.fromImmuneNetwork(immuneNetwork)

        session.add(network)
        session.commit()

        session.refresh(network)

        return network.network_id
        
def get_network(network_id):

    with Session(engine) as session:

        return session.query(Network).options(
            selectinload(Network.network_edges),
            selectinload(Network.source_repertoire)
            .selectinload(Repertoire.clonotypes)
        ).filter(
            Network.network_id == network_id
        ).first()



