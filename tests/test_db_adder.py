from pvfsga.db import adder
from pvfsga.db.orm.molecule import Molecule as DBChem, Scaffold, RGroup
from pvfsga.db.orm.reaction import (
    Reaction,
    ReactionStep,
    ReactionProcedure as DBReactionProcedure,
)
from pvfsga.models import ReactionProcedure


def test_add_reaction_procedure(session):
    """Tests that when a reaction procedure is added it also adds the steps and
    reactions and doesn't add the procedure if it already exists
    """
    assert (
        session.query(DBReactionProcedure).filter_by(name="TestReaction").count() == 0
    )
    procedure = ReactionProcedure(
        name="TestReaction",
        steps=[
            "[C:6]1[C:1][O:2][C:3](=[O:4])[C:5]1.[S:7]>>[C:6]1[C:1][O:2][C:3](=[S:7])[C:5]1.[O:4]",
            "[C:6]1[C:1][O:2][C:3](=[S:7])[C:5]1.[O:4]>>[*][C:3](=[O:2])[C:5][C:6][C:1][S:7][*].[O:4]",
        ],
    )
    procedure_id = adder.add_reaction_procedure_and_or_get_id(procedure, session)
    assert (
        session.query(DBReactionProcedure).filter_by(name="TestReaction").count() == 1
    )
    # confirm it isn't readded
    adder.add_reaction_procedure_and_or_get_id(procedure, session)
    assert (
        session.query(DBReactionProcedure).filter_by(name="TestReaction").count() == 1
    )
    assert (
        session.query(ReactionStep)
        .filter_by(reaction_procedure_id=procedure_id)
        .count()
        == 2
    )
    assert (
        session.query(Reaction)
        .filter_by(
            smarts="[C:6]1[C:1][O:2][C:3](=[O:4])[C:5]1.[S:7]>>[C:6]1[C:1][O:2][C:3](=[S:7])[C:5]1.[O:4]"
        )
        .count()
        == 1
    )
    assert (
        session.query(Reaction)
        .filter_by(
            smarts="[C:6]1[C:1][O:2][C:3](=[S:7])[C:5]1.[O:4]>>[*][C:3](=[O:2])[C:5][C:6][C:1][S:7][*].[O:4]"
        )
        .count()
        == 1
    )


def test_add_procedure_with_same_reaction(session):
    """ "Tests that when a reaction procedure that has the same reaction as another
    is added, both use the same reaction id in the db"""
    procedure1 = ReactionProcedure(
        name="TestReaction",
        steps=[
            "[C:6]1[C:1][O:2][C:3](=[O:4])[C:5]1.[S:7]>>[C:6]1[C:1][O:2][C:3](=[S:7])[C:5]1.[O:4]",
            "[C:6]1[C:1][O:2][C:3](=[S:7])[C:5]1.[O:4]>>[*][C:3](=[O:2])[C:5][C:6][C:1][S:7][*].[O:4]",
        ],
    )
    procedure2 = ReactionProcedure(
        name="TestReaction2",
        steps=[
            "[C:6]1[C:1][O:2][C:3](=[O:4])[C:5]1.[S:7]>>[C:6]1[C:1][O:2][C:3](=[S:7])[C:5]1.[O:4]"
        ],
    )

    procedure1_id = adder.add_reaction_procedure_and_or_get_id(procedure1, session)
    procedure2_id = adder.add_reaction_procedure_and_or_get_id(procedure2, session)
    assert (
        session.query(Reaction)
        .filter_by(
            smarts="[C:6]1[C:1][O:2][C:3](=[O:4])[C:5]1.[S:7]>>[C:6]1[C:1][O:2][C:3](=[S:7])[C:5]1.[O:4]"
        )
        .count()
        == 1
    )
    assert (
        session.query(Reaction)
        .filter_by(
            smarts="[C:6]1[C:1][O:2][C:3](=[S:7])[C:5]1.[O:4]>>[*][C:3](=[O:2])[C:5][C:6][C:1][S:7][*].[O:4]"
        )
        .count()
        == 1
    )
    assert session.query(Reaction).count() == 2
    assert (
        session.query(DBReactionProcedure)
        .filter(
            DBReactionProcedure.reaction_procedure_id.in_(
                (procedure1_id, procedure2_id)
            )
        )
        .count()
        == 2
    )


def test_add_rgroup(session):
    rgroup = "CCOCC"
    adder.add_rgroup_and_or_get_id(rgroup, session)
    assert session.query(RGroup).filter_by(smiles=rgroup).count() == 1
    # confirm rgroup isn't readded
    adder.add_rgroup_and_or_get_id(rgroup, session)
    assert session.query(RGroup).filter_by(smiles=rgroup).count() == 1


def test_add_scaffold(session):
    scaffold = "O=C1CC([*:1])CO([*:2])1"
    adder.add_scaffold_and_or_get_id(scaffold, session)
    assert session.query(Scaffold).filter_by(structure=scaffold).count() == 1
    # confirm scaffold isn't readded
    adder.add_scaffold_and_or_get_id(scaffold, session)
    assert session.query(Scaffold).filter_by(structure=scaffold).count() == 1


def test_add_chem(session):
    smiles = "CCOCC"
    scaffold_id = 10000
    adder.add_molecule_and_get_id(smiles, scaffold_id, session)
    assert session.query(DBChem).filter_by(smiles=smiles).count() == 1
    # confirm that we aren't checking to prevent degeneracy here
    adder.add_molecule_and_get_id(smiles, scaffold_id, session)
    assert session.query(DBChem).filter_by(smiles=smiles).count() == 2
