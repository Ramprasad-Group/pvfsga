from typing import List

import pytest

from pvfsga.models import Polymer
from .test_crossover import (
    obj_2,  # noqa: F401
    silly_chem_1,  # noqa: F401
    silly_chem_2,  # noqa: F401
    silly_chem_3,  # noqa: F401
)


def fp(polymers: List[Polymer]) -> List[Polymer]:
    for polymer in polymers:
        polymer.fingerprint = {"fp1": "test"}
    return polymers


@pytest.fixture
def list_molecules(silly_chem_1, silly_chem_2):  # noqa: F811
    return [silly_chem_1, silly_chem_2]


def test_exchange(obj_2, silly_chem_1):  # noqa: F811
    new_molecule = obj_2.exchange(silly_chem_1, 3, "[*]COC")
    assert new_molecule.r_groups["R3"] == "[*]COC"
    new_molecule = obj_2.exchange(silly_chem_1, 1, "")
    assert new_molecule.r_groups["R1"] == ""
    new_molecule = obj_2.exchange(silly_chem_1, 3, "[*]COC")
    assert new_molecule.r_groups["R3"] == "[*]COC"


def test_mutate(obj_2, silly_chem_2):  # noqa: F811
    silly_chem_2_old = silly_chem_2.copy()
    flag = False
    for _ in range(10):
        mutated_molecule = obj_2.mutate(silly_chem_2)
        if mutated_molecule != silly_chem_2_old:
            flag = True
    assert flag


def test_mutate_all(obj_2, list_molecules):  # noqa: F811
    list_molecules_old = [mol.copy() for mol in list_molecules]
    new_mutated_list = obj_2.mutate_all(list_molecules, 1)
    diff = False
    for i in range(len(list_molecules_old)):
        m1 = list_molecules_old[i]
        m2 = new_mutated_list[i]
        if m1.smiles != m2.smiles:
            diff = True
            break
    assert diff
