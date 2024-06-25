from rdkit.Chem import rdChemReactions
from rdkit import Chem
import pandas as pd

from pvfsga.db.orm.polymer import Polymer as Polymer_DB
from pvfsga.db.orm.polymer import PolymerProperty, PolymerFingerprint
from pvfsga.db.orm.molecule import Molecule as Molecule_DB
from pvfsga.db.orm.molecule import Scaffold, RGroup, RGroupMapping
from pvfsga.db.orm.reaction import Reaction, ReactionStep, ReactionProcedure


def test_get_pol(session):
    """Simple test to see if connection established and that adding datetime works"""
    pol = session.query(Polymer_DB).first()
    assert pol.smiles == "*SCC(CC(*)=O)N=C=S"


def test_get_pol_property(session):
    prop = session.query(PolymerProperty).filter(PolymerProperty.value > 200).all()
    assert len(prop) == 1
    assert prop[0].value == 350
    assert prop[0].prop == "Glass Transition Temperature"


def test_get_pol_fingerprint(session):
    """Tests if fingerprints properly loaded and tests converting from long table into
    wide table
    """
    fingerprints = (
        session.query(
            PolymerFingerprint.key, PolymerFingerprint.value, Polymer_DB.smiles
        )
        .join(PolymerFingerprint)
        .all()
    )
    # Convert fingerprints from long table to wide table
    data = {}
    for key, value, polymer in fingerprints:
        if polymer not in data:
            data[polymer] = {"polymer": polymer}
        data[polymer][key] = value
    df = pd.DataFrame(data.values())
    for col in ["fp1", "fp2", "fp3", "polymer"]:
        assert col in df.columns


def test_get_molecule_scaffold_and_rgroup(session):
    chem, scaffold, rgroup = (
        session.query(Molecule_DB, Scaffold, RGroup)
        .join(Molecule_DB.scaffold)
        .join(Molecule_DB.rgroup_mappings)
        .join(RGroupMapping.rgroup)
        .filter(Molecule_DB.molecule_id == 10000)
        .first()
    )
    assert chem.smiles == "O=C1CC(N=C=S)CO1"
    assert scaffold.structure == "O=C1CC([*:1])CO1"
    assert rgroup.smiles == "N=C=S"


def test_get_reaction_procedure(session):
    steps = (
        session.query(Reaction.smarts, ReactionStep.step)
        .join(ReactionStep.reaction)
        .join(ReactionStep.procedure)
        .filter(ReactionProcedure.reaction_procedure_id == 10000)
        .order_by(ReactionStep.step)
        .all()
    )
    assert steps[0][0] == (
        "[C:6]1[C:1][O:2][C:3](=[O:4])[C:5]1.[S:7]"
        + ">>[C:6]1[C:1][O:2][C:3](=[S:7])[C:5]1.[O:4]"
    )
    assert steps[1][0] == (
        "[C:6]1[C:1][O:2][C:3](=[S:7])[C:5]1.[O:4]>>"
        + "[*][C:3](=[O:2])[C:5][C:6][C:1][S:7][*].[O:4]"
    )


def test_run_reaction(session):
    chem = (
        session.query(Molecule_DB.smiles)
        .filter(Molecule_DB.molecule_id == 10000)
        .first()[0]
    )
    steps = (
        session.query(Reaction.smarts, ReactionStep.step)
        .join(ReactionStep.reaction)
        .join(ReactionStep.procedure)
        .filter(ReactionProcedure.reaction_procedure_id == 10000)
        .order_by(ReactionStep.step)
        .all()
    )
    products = (
        Chem.MolFromSmiles(chem),
        Chem.MolFromSmiles("S"),
    )
    for step in steps:
        smarts = step[0]
        rxn = rdChemReactions.ReactionFromSmarts(smarts)
        products = rxn.RunReactants(products)[0]
    polymer = Chem.MolToSmiles(products[0])
    assert session.query(Polymer_DB).filter(Polymer_DB.smiles == polymer).count() > 0
