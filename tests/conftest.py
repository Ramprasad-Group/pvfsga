from pathlib import Path
import json

import pytest
from sqlalchemy_utils import database_exists, create_database, drop_database
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session

from pvfsga.db import Base
from pvfsga.db.orm.polymer import Polymer as Polymer_DB
from pvfsga.db.orm.polymer import (
    PolymerFingerprint,
    PolymerProperty,
    PolymerGeneration,
)
from pvfsga.db.orm.molecule import Molecule as Molecule_DB

from pvfsga.db.orm.molecule import Scaffold, RGroup, RGroupMapping
from pvfsga.db.orm.reaction import (
    Reaction,
    ReactionStep,
    ReactionProcedure,
    ReactionPolymerMapping,
    ReactionScaffoldMapping,
)


@pytest.fixture(scope="session")
def engine():
    db_server = "sqlite:///test_run.sqlite"
    engine = create_engine(db_server)
    if not database_exists(engine.url):
        create_database(engine.url)
    else:
        raise FileExistsError(
            "Database already exists. Since the database is "
            + "dropped at the end of testing, stopping the test."
        )
    yield engine

    drop_database(engine.url)


@pytest.fixture(scope="module")
def setup_database(engine):
    Base.metadata.create_all(engine)

    yield

    Base.metadata.drop_all(engine)


@pytest.fixture(scope="function")
def session(engine, setup_database):
    connection = engine.connect()
    transaction = connection.begin()
    session = scoped_session(
        sessionmaker(autocommit=False, autoflush=False, bind=connection)
    )
    try:
        seed_database(session)
    except Exception as e:
        print(e)
        transaction.rollback()
        return
    yield session
    session.close()
    transaction.rollback()
    connection.close()


def seed_database(session):
    """Seeds database for testing"""
    folder = Path(__file__).parent.absolute() / "test_db_seeds"

    with open(str(folder / "polymer.json"), "r") as fin:
        module = json.loads(fin.read())  # module: Dict
    for data in module["Polymer"]:
        session.add(Polymer_DB(**data))
        session.commit()
    for data in module["PolymerFingerprint"]:
        session.add(PolymerFingerprint(**data))
        session.commit()
    for data in module["PolymerProperty"]:
        session.add(PolymerProperty(**data))
        session.commit()
    for data in module["PolymerGeneration"]:
        session.add(PolymerGeneration(**data))
        session.commit()

    with open(str(folder / "molecule.json"), "r") as fin:
        module = json.loads(fin.read())
    for data in module["Molecule"]:
        session.add(Molecule_DB(**data))
        session.commit()
    for data in module["Scaffold"]:
        session.add(Scaffold(**data))
        session.commit()
    for data in module["RGroup"]:
        session.add(RGroup(**data))
        session.commit()
    for data in module["RGroupMapping"]:
        session.add(RGroupMapping(**data))
        session.commit()

    with open(str(folder / "reaction.json"), "r") as fin:
        module = json.loads(fin.read())
    for data in module["Reaction"]:
        session.add(Reaction(**data))
        session.commit()
    for data in module["ReactionStep"]:
        session.add(ReactionStep(**data))
        session.commit()
    for data in module["ReactionProcedure"]:
        session.add(ReactionProcedure(**data))
        session.commit()
    for data in module["ReactionPolymerMapping"]:
        session.add(ReactionPolymerMapping(**data))
        session.commit()
    for data in module["ReactionScaffoldMapping"]:
        session.add(ReactionScaffoldMapping(**data))
        session.commit()
