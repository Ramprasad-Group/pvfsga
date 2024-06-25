"""Python script that runs the GA, and includes the fitness functions"""
import logging
from pvfsga import store_poly_chem
from pvfsga.models import Polymer
import numpy as np
from pvfsga.db import adder
from sklearn.preprocessing import MinMaxScaler


logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.propagate = False


def run(
    obj_GA,
    predictor: callable,
    selected_species_number: int = 50,
    fitness_function: str = "default",
    props_to_fit: dict[str:dict] = {
        "Tg": {"target": 342.0, "max": 600.00, "min": 373, "cap": True}
    },
    generation_threshold: int = 10,
):
    """Calls the recursive function to run the genetic algorithm, and then connects
    to the output file to process the data

    Attributes:
        obj_GA(pvfsga.main.GA_society):
            An instance of the the GA class that has the reproduction function

        predictor (callable):
            A callable object representing the predictor function, that takes in a
            list[Polymers] and returns list[Polymers] with Polymer.properties

        selected_species_number (int, optional):
            The number of selected species. Defaults to 50.

        fitness_function (str):
            string identifying which fitness function to utilize. Defualt: default

        props_to_fit (dict[str:dict]):
            dictionary containing the desired properties to optimize,
            along with the target values


    Returns:
        None
    """
    parent_fitness_results = []
    parent_0 = [
        pol
        for pol in obj_GA.polymerize_all(obj_GA.initial_population.copy())
        if pol is not None
    ]
    parent_0 = obj_GA.fingerprint_efficient(parent_0)
    store_poly_chem.store_polymer(obj_GA, parent_0)
    if fitness_function == "default":
        parent_fitness_results = predictor(parent_0, props_to_fit.keys())
        add_predicted_to_db(obj_GA, parent_fitness_results)

        recursive_func(
            obj_GA,
            predictor,
            parent=parent_fitness_results,
            selected_species_number=selected_species_number,
            fitness_function=fitness_function,
            props_to_fit=props_to_fit,
            generation_threshold=generation_threshold,
        )


def recursive_func(
    obj_GA,
    predictor: callable,
    parent: list[Polymer] = None,
    selected_species_number: int = 50,
    fitness_function: str = "default",
    props_to_fit: dict[str:dict] = {
        "Tg": {"target": 342.0, "max": 600.00, "min": 373, "cap": True}
    },
    generation_threshold: float = 0.9,
) -> list[Polymer]:
    """The recursive function used to create generations and pick the top
    species generation based on a certain fitness function

    Attributes:
        parent (list[Polymers]):
            The list of all parents (fitted and sorted) at a specific generation.


    Returns:
        fitted_children (list[Polymer]):
            a list of the species desired after reaching a
            specific fitness value
    """
    if obj_GA.verbose > 0:
        print("This is generation number", obj_GA.generation)
    children = obj_GA.reproduce(parent)
    polymer_list = obj_GA.fingerprint_efficient(children)
    store_poly_chem.store_polymer(obj_GA, polymer_list)

    if obj_GA.verbose > 0:
        print("this is children num", len(children))
    if fitness_function == "default":
        children_fitness_results = fitness_evalutaion(predictor, children, props_to_fit)
        add_predicted_to_db(obj_GA, children_fitness_results)
        sorted_children_fitness_results = sorted(
            children_fitness_results,
            key=lambda x: x.fitness["default"],
            reverse=True,
        )

        # The Update happens here
        fitted_children = sorted_children_fitness_results[:selected_species_number]
        values = [d.fitness["default"] for d in fitted_children]
        average = sum(values) / len(values) if len(values) > 0 else 0
        print("this is avg", average)
        if obj_GA.generation > generation_threshold:
            return fitted_children
        else:
            return recursive_func(
                obj_GA,
                predictor,
                parent=fitted_children,
                selected_species_number=selected_species_number,
                fitness_function=fitness_function,
                props_to_fit=props_to_fit,
                generation_threshold=generation_threshold,
            )

    else:
        raise ValueError("Please choose a valid fitness function")


def fitness_evalutaion(
    predictor: callable,
    children: list[Polymer],
    props_to_fit: dict[str:dict],
):
    """evalutes the fitness based on the criteria listed in props_to_fit

    Attributes:
        children (list[Polymers]):
            The list of all children to have the fitness evaluated

    Returns:
        children_fitness_results (list[Polymer]):
            a list of the polymers with their fitness evaluation
    """

    children_fitness_results = predictor(children, list(props_to_fit.keys()))
    scaled_dict = {}
    for prop in props_to_fit.keys():
        identical_flag = False

        vals = np.array([poly.properties[prop] for poly in children_fitness_results])

        if props_to_fit[prop]["cap"]:
            vals = np.array(
                [max(value, props_to_fit[prop]["target"]) for value in vals]
            )
            vals = -1 * vals
        else:
            vals = np.array(
                [min(value, props_to_fit[prop]["target"]) for value in vals]
            )

        if np.all(np.isclose(vals, vals[0], atol=1e-8)):
            identical_flag = True

        vals = np.reshape(vals, (len(vals), -1))
        scaler = MinMaxScaler()
        scaled_vals = scaler.fit_transform(vals)
        scaled_vals = np.array(scaled_vals).flatten()
        if identical_flag:
            logger.info(f"All values in {prop} have been optimized")
            print("this is flag")
            scaled_vals = 1 - scaled_vals

        scaled_dict[prop] = scaled_vals

    for i, poly in enumerate(children_fitness_results):
        fit = 0
        for prop in list(props_to_fit.keys()):
            fit += 1 / (len(list(props_to_fit.keys()))) * scaled_dict[prop][i]
        fit_dict = {"default": fit}
        poly.fitness = fit_dict
        poly.properties["fitness"] = fit

    return children_fitness_results


def add_predicted_to_db(obj_GA, polymer_list: list[Polymer]):
    """Adds predicted polymer properties to the database.

    Attributes:
        obj_GA(pvfsga.main.GA_society):
            An instance of the the GA class that has the reproduction function

        polymer_list (list[Polymer]):
            A list of Polymer objects with predicted properties to be added to the DB

    Returns:
        None
    """
    polymer_dict_temp = set()
    with obj_GA.Session() as session:
        for iterated_polymer in polymer_list:
            if iterated_polymer.canonical_smiles not in obj_GA.polymer_property:
                for dummy_prop, dummy_prop_value in iterated_polymer.properties.items():
                    adder.add_polymer_property(
                        pol_id=obj_GA.polymer_dict[iterated_polymer.canonical_smiles],
                        prop=dummy_prop,
                        prop_unit="unit",
                        value=dummy_prop_value,
                        session=session,
                    )
                    obj_GA.polymer_property[
                        iterated_polymer.canonical_smiles
                    ] = iterated_polymer.properties
                    polymer_dict_temp.add(
                        obj_GA.polymer_dict[iterated_polymer.canonical_smiles]
                    )
                    session.commit()
