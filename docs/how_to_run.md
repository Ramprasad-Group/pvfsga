## Installation for Linux / macOS / Windows
### [Poetry](https://python-poetry.org/)

```bash
poetry add git+ssh://git@github.com/jdkern11/pvfsga.git#main
```

To specify a specific release, change main to the version number (v#.#.# e.g., v1.1.1). 

Example:
```bash
poetry add git+ssh://git@github.com/jdkern11/pvfsga.git#v1.1.0
```
Once how you have setup this package with poetry, you can then use the below script to run the GA within the the poetry environment just established. 
## How to run

```Python
"""Python file used to manage inputs, fitness functions, and more"""
import joblib
import logging
import sys
#import list


from pvfsga import main, run_model
from pvfsga.models import Polymer

import pandas as pd
from rdkit import RDLogger, Chem
from rdkit.Chem import AllChem
from rdkit.Contrib.SA_Score.sascorer import calculateScore



RDLogger.DisableLog("rdApp.*")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s",
    filename="logfile.log",
    filemode="w",
)
logging.getLogger("pgfingerprinting.fp").setLevel(logging.ERROR)


logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.propagate = False

def fingerprint_all(polymers: list[Polymer]):
    """
    Fingerprint all polymers in the given list.

    Args:
        polymers (list[Polymer]): List of Polymer objects.

    Returns:
        list[Polymer]: List of polymers with fingerprint attribute.
    """
    prog = 0

    pols_to_return = []
    for polymer in polymers:
        PROG = 100 * round(prog / (len(polymers)), 3)
        sys.stdout.write("fingerprinting \r%d%% done" % PROG)
        sys.stdout.flush()
        prog += 1

        m1 = Chem.MolFromSmiles(polymer.smiles)
        bi = {}
        fpgen = AllChem.GetMorganFingerprintAsBitVect(m1,radius=2,bitInfo=bi, nBits=1024)
        bit_list = list(fpgen)
        fingerpint_dict= {str(index): value for index, value in enumerate(bit_list)}
        polymer.fingerprint = fingerpint_dict.copy()
        #print( polymer.fingerprint)
        pols_to_return.append(polymer)
        
    return pols_to_return

        

def predict_all(polymers: list[Polymer], prop_to_find):
    """
    Predict properties for all polymers in the given list. 
    Transform property values as necessary to minimize (e.g. SA Score),
     or attain a value within a range (e.g. Enthalpy)

    Args:
        polymers (list[Polymer]): List of Polymer objects.

    Returns:
        list[Polymer]: List of polymers with properties attribute.
    """
    dh_reg = joblib.load('fakedh_linear_model.pkl')
    Tg_reg = joblib.load('fakeTg_linear_model.pkl')
    polymers_predicted = []
    prog = 0

    for polymer in polymers:
        PROG = 100 * round(prog / (len(polymers)), 3)
        sys.stdout.write("pred \r%d%% done" % PROG)
        sys.stdout.flush()
        prog += 1
        properties_temp  = {}
        properties_temp["sa_score"] = (-1)*calculateScore(Chem.MolFromSmiles(polymer.molecule.smiles))
        properties_temp["fake_Tg"] = Tg_reg.predict([list(polymer.fingerprint.values())])[0]
        properties_temp["fake_dH"] = dh_reg.predict([list(polymer.fingerprint.values())])[0]
        enthalpy_val =  properties_temp["fake_dH"]
        if enthalpy_val < -20:
            properties_temp["fake_dH_eff"] = (
               enthalpy_val + 30
            )
        elif enthalpy_val >= -20 and enthalpy_val <= -10:
           properties_temp["fake_dH_eff"] = 10
        elif enthalpy_val >= -10:
            properties_temp["fake_dH_eff"] = (
                -1 * enthalpy_val
            )
        polymer.properties  = properties_temp.copy()
        polymers_predicted.append(polymer)
    

    return polymers_predicted


######### R-groups file ##########################
csv_file_path = "small_rgroups.csv"
df = pd.read_csv(csv_file_path)
df["rgroups"].fillna("", inplace=True)

# Convert the "R-groups" column of the DataFrame to a list
r_groups = list(set(df["rgroups"].tolist()))
print(len(r_groups))
#########Reactions file ##########################


reaction_procedures = {
    "ROMPs": {
        "ROMP ROP": [
            "[C;R:1]=[C;R:2].[Ge:3]>>[C;R:1]=[*][Ge:3][C;R:2]",
            "[Ge:3][C;R:2]>>[*]=[C:2].[Ge:3]",
        ]
    },
}

reaction_dict = reaction_procedures["ROMPs"]
######### Scaffolds ##########################
scaffolds = [" [*:1]C1C2C=CC(C2)C1[*:2]","[*:1]C1C=CC([*:2])C1","[*:1]C1([*:2])CC2C=CC1C2"] #List containing all scaffolds for a single run


scaffold_reaction_dict = {}
for scaffold in scaffolds:
    scaffold_reaction_dict[scaffold] = list(reaction_dict.keys())
    

############### User-Defined Properties############
properties_target = {
    "sa_score": {"target": -3, "cap": False},
    "fake_Tg": {"target": 373, "cap": False},
    "fake_dH_eff": {"target": 10, "cap": False},
} # Dictionary that contains the property name, target, and whether the property should be above or below the target value



obj = main.GA_society(
    r_groups,
    reaction_dict,
    scaffold_reaction_dict,
    initial_population_num=100,
    num_families=50,
    num_children_per_family=3,
    partner_selection="random",
    fp_function=fingerprint_all,
    verbose = 1
)

try:
    run_model.run(
            obj,
            predict_all,
            selected_species_number=100,
            fitness_function="default",
            props_to_fit=properties_target,
            generation_threshold=100,
        )
except Exception as e:
    logger.error(f"GA run failed with error: {e}")
```

1. initial_population_num: The number of molecules/polymers generated for the first population
2. num_families: the number of families (couples) the polymers form each generation. This couple then generates other monomers/polymers through cross-over and mutations. 
3. num_children_per_family: the number of monomer children each couple creates when they are paired up
4. partner_selection: the method used to choose partners. Currently, only the 'Random' scheme is implemented. This means that, from the top polymers chosen each generation, a certain number of couples, specified by 'num_families,' are selected randomly to form a subset of polymers.
5. selected_species_number: The number of top performing polymers chosen each generation to be mated together. These polymers are paired up to form "num_families" couples. 
6. generation_threshold: The number of generations the GA will run.
