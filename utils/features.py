import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Lipinski
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from rdkit.Chem import AllChem
import combichem.utils.io as io


def get_filter_catalog(mol):
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.ZINC)

    filters = FilterCatalog.FilterCatalog(params)
    if filters.HasMatch(mol):
        decription = []
        reference = []
        scope = []

        matches = filters.GetMatches(mol)

        for match in matches:
            decription.append(match.GetDescription())
            reference.append(match.GetProp("Reference"))
            scope.append(match.GetProp("Scope"))
        return {"Descriptors":"; ".join(decription), "Reference":"; ".join(reference), "Score":"; ".join(scope)}

    else:
        return {"Descriptors":"", "Reference":"", "Scope":""}


def mmff_energy(mol):
    try:
        mp = AllChem.MMFFGetMoleculeProperties(mol)
        ffm = AllChem.MMFFGetMoleculeForceField(mol, mp)
        energy_value_M = ffm.CalcEnergy()
        return energy_value_M
    except:
        return "None"

def get_descriptors(mol):
    props_names = ["exactmw",
                  "NumRotatableBonds",
                  "NumHBD",
                  "NumHBA",
                  "NumRings",
                  "CrippenClogP"]

    props = rdMolDescriptors.Properties(props_names)
    props = zip(props.GetPropertyNames(),
                props.ComputeProperties(mol))

    props = {name:value for name, value in props}
    props["HeavyAtoms"] = Lipinski.HeavyAtomCount(mol)

    return props

def standardization(mol_original):
    mol = Chem.Mol(mol_original)

    Chem.SanitizeMol(mol)
    Chem.rdmolops.SetAromaticity(mol, model=Chem.rdmolops.AromaticityModel.AROMATICITY_MDL)
    Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)

    smiles = io.write_smi(mol, remove_hs=False)
    descriptors = get_descriptors(mol)
    filters = get_filter_catalog(mol)
    energy = mmff_energy(mol)

    dic = {"smi": smiles, "MMFF": energy}

    for desc, value in descriptors.items():
        dic[desc] = value

    for filt, value in filters.items():
        dic[filt] = value

    return dic


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", type=str, required=True)
    parser.add_argument("-header", type=bool, default=True)
    
    args = parser.parse_args()
    input_mol = args.input
    header = args.header

    mol = io.read_sdf(input_mol, add_hs=False, remove_hs=False)[0]
    props = standardization(mol)

    if header:
        terms = ['smi', 'MMFF', 'exactmw', 'NumRotatableBonds', 'NumHBD', 'NumHBA', 'NumRings', 'CrippenClogP', 'HeavyAtoms', 'Descriptors', 'Reference', 'Scope']

        print("Input\t" + "\t".join(terms))

    vals = [input_mol] + list(props.values())
    vals = "\t".join([str(elem) for elem in vals])

    print(vals)
    

if __name__ == "__main__":
    main()