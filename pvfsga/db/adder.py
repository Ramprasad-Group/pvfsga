from typing import Dict, Any
import logging

from pvfsga.db.orm.polymer import (
    Polymer,
    PolymerProperty,
    PolymerFingerprint,
    PolymerGeneration,
)
from pvfsga.db.orm.molecule import Molecule, Scaffold, RGroup, RGroupMapping
from pvfsga.db.orm.reaction import (
    Reaction,
    ReactionStep,
    ReactionProcedure as DBReactionProcedure,
    ReactionPolymerMapping,
    ReactionScaffoldMapping,
)
from pvfsga.models import ReactionProcedure

logger = logging.getLogger(__name__)


def add_reaction_scaffold_mappings(scaffold_id, reaction_procedure_id, session):
    dat = {"scaffold_id": scaffold_id, "reaction_procedure_id": reaction_procedure_id}
    add_row_and_get_obj(dat, ReactionScaffoldMapping, session)


def add_polymer_property(pol_id, prop, prop_unit, value, session):
    dat = {"pol_id": pol_id, "prop": prop, "prop_unit": prop_unit, "value": value}
    add_row_and_get_obj(dat, PolymerProperty, session)


def add_reaction_polymer_mapping(reaction_procedure_id, molecule_id, pol_id, session):
    dat = {
        "reaction_procedure_id": reaction_procedure_id,
        "molecule_id": molecule_id,
        "pol_id": pol_id,
    }
    add_row_and_get_obj(dat, ReactionPolymerMapping, session)


def add_polymer_fingerprint(key, value, pol_id, session):
    dat = {"key": key, "value": value, "pol_id": pol_id}
    add_row_and_get_obj(dat, PolymerFingerprint, session)


def add_polymer_generation(gen, pol_id, session):
    dat = {"generation": gen, "pol_id": pol_id}
    add_row_and_get_obj(dat, PolymerGeneration, session)


def add_polymer_and_get_id(smiles, canonical_smiles, session) -> int:
    dat = {"smiles": smiles, "canonical_smiles": canonical_smiles}
    return add_row_and_get_obj(dat, Polymer, session).pol_id


def add_rgroup_mapping(molecule_id, rgroup_id, index, session):
    dat = {"molecule_id": molecule_id, "rgroup_id": rgroup_id, "index": index}
    add_row_and_get_obj(dat, RGroupMapping, session)


def add_molecule_and_get_id(smiles: str, scaffold_id: int, session) -> int:
    # Not checking if molecule already added here, as we expect it to be slow
    # this means check has to be done prior to calling this function
    dat = {"smiles": smiles, "scaffold_id": scaffold_id}
    return add_row_and_get_obj(dat, Molecule, session).molecule_id


def add_reaction_procedure_and_or_get_id(procedure: ReactionProcedure, session) -> int:
    """Adds reaction procedure and returns its reaction_procedure_id"""
    dat = {"name": procedure.name}

    query = session.query(DBReactionProcedure).filter_by(**dat)
    if query.count() > 0:
        logger.info(f"{procedure} already added")
        return query.first().reaction_procedure_id
    else:
        proc = DBReactionProcedure(**dat)
        session.add(proc)
        session.flush()
        session.refresh(proc)
        for index, step in enumerate(procedure.steps):
            reaction_id = add_reaction_and_or_get_id(step, session)
            add_reaction_step(
                proc.reaction_procedure_id, reaction_id, index + 1, session
            )
        return proc.reaction_procedure_id


def add_reaction_and_or_get_id(smarts: str, session) -> int:
    """Adds reaction to database and returns its id"""
    dat = {"smarts": smarts}
    return add_row_if_count_zero_and_get_obj(dat, Reaction, session).reaction_id


def add_reaction_step(reaction_procedure_id: int, reaction_id: int, step: int, session):
    dat = {
        "reaction_procedure_id": reaction_procedure_id,
        "reaction_id": reaction_id,
        "step": step,
    }
    add_row_if_count_zero_and_get_obj(dat, ReactionStep, session)


def add_rgroup_and_or_get_id(smiles: str, session) -> int:
    dat = {"smiles": smiles}
    return add_row_if_count_zero_and_get_obj(dat, RGroup, session).rgroup_id


def add_scaffold_and_or_get_id(structure: str, session) -> int:
    dat = {"structure": structure}
    return add_row_if_count_zero_and_get_obj(dat, Scaffold, session).scaffold_id


def add_row_if_count_zero_and_get_obj(dat: Dict[str, Any], Table, session) -> Any:
    """Adds row of data to general table object from db orm and returns the added
    and flushed row
    """
    query = session.query(Table).filter_by(**dat)
    if query.count() > 0:
        logger.info(f"{Table} {dat} already added")
        return query.first()
    row = Table(**dat)
    session.add(row)
    session.flush()
    session.refresh(row)
    return row


def add_row_and_get_obj(dat: Dict[str, Any], Table, session) -> Any:
    """Adds row of data to general table object from db orm and returns the added
    and flushed row
    """
    row = Table(**dat)
    session.add(row)
    session.flush()
    session.refresh(row)
    return row
