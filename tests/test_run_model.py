import pytest

from pvfsga import run_model
from pvfsga.models import Polymer, ReactionProcedure
from .test_crossover import obj_2  # noqa: F401


def fingerprint_all(polymers: list[Polymer]):
    """
    Fingerprint all polymers in the given list.

    Args:
        polymers (list[Polymer]): List of Polymer objects.

    Returns:
        list[Polymer]: List of polymers with fingerprint attribute.
    """
    for polymer in polymers:
        polymer.fingerprint = {"test_fp": 1}
    return polymers


def predict_all(polymers: list[Polymer], prop_to_fit):
    """
    Predict properties for all polymers in the given list.

    Args:
        polymers (list[Polymer]): List of Polymer objects.

    Returns:
        list[Polymer]: List of polymers with properties attribute.
    """
    for polymer in polymers:
        if polymer.smiles is not None:
            polymer.properties = {"length": len(polymer.smiles)}
        else:
            polymer.properties = {"length": len(polymer.smiles)}
    return polymers


@pytest.fixture
def reaction_procedure():
    return ReactionProcedure(
        name="test1",
        steps=[
            "[C:1]1[O:2][C:3](=[S:4])[C:5][C:6]1>>[*][S:4][C:1][C:6][C:5][C:3](=[O:2])[*]"
        ],
    )


def test_recursive_func(obj_2):  # noqa: F811
    parent_0 = obj_2.polymerize_all(obj_2.initial_population)
    fitted = run_model.recursive_func(
        obj_2,
        predict_all,
        parent_0,
        selected_species_number=3,
        fitness_function="default",
        props_to_fit={"length": {"target": 5, "max": 25, "min": 0, "cap": True}},
        generation_threshold=1,
    )
    sum = 0

    for poly in fitted:
        sum = sum + poly.properties["length"]
    average = sum / len(fitted)
    assert isinstance(fitted[0], Polymer)
    assert average > 0
