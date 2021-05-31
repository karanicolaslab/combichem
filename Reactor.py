from itertools import product
from rdkit import Chem
from rdkit.Chem import AllChem
import utils.io as io

def run_reactor(reactants, reaction, add_hs=True):
    """Perform reaction with reactants

    Args:
        reactant (list): List of reactants
        reaction (str): SMARTS reaction

    Returns:
        set: List of products of reaction
    """
    rxn = AllChem.ReactionFromSmarts(reaction)

    if rxn is None:
        raise Exception("Cannot read {0}".format(reaction))

    for i, elem in enumerate(reactants):
        pattern = reaction.split(">>")[0].split(".")[i]
        pattern = Chem.MolFromSmarts(pattern)
        if elem.HasSubstructMatch(pattern) is False:
            print(Chem.MolToSmiles(elem), Chem.MolToSmiles(pattern))
            raise Exception("Reactant {0} doesn't match to template".format(i+1))

    product = rxn.RunReactants(reactants)

    prods = []

    for i in range(len(product)):
        mol = product[i][0]
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        if add_hs:
            mol = io.add_hydrogens(mol)
        prods.append(Chem.MolToSmiles(mol))

    return list(set(prods))



def reactor(reactants, reaction, prefix="R1"):
    
    
    bbs = []

    for i, file in enumerate(reactants):
        with open(file, "r") as f:
            f = f.read().split("\n")
            f = [elem for elem in f if elem != ""]
        bbs.append(f)

    # parse reaction
    with open(reaction, "r") as f:
        reaction = f.read().split("\n")[0]

    bbs = list(product(*bbs))

    products = []

    # perform reaction
    for i, line in enumerate(bbs):
        smi_mols = []
        smi_names = []
        for i, smi in enumerate(line):
            if "\t" in smi:
                smi_mol, smi_name = smi.split("\t")
            elif "," in smi:
                smi_mol, smi_name = smi.split(",")
            elif " " in smi:
                smi_mol, smi_name = smi.split(" ")
            else:
                smi_mol = smi
                smi_name = "Mol_" + str(i+1)

            smi_mols.append(smi_mol)
            smi_names.append(smi_name)

        smi_mols_new = []

        try:
            for smi in smi_mols:
                smi_mols_new.append(io.read_smi(smi, add_hs=False))
        except:
            continue

        smi_mols = smi_mols_new[:]

        try:
            prods = run_reactor(smi_mols, reaction, add_hs=False)

            for j, prod in enumerate(prods):
                products.append([prod, prefix + "_" + "_".join(smi_names) + "v" + str(j)])

        except Exception as e:
            print(e)
            # print(io.write_smi(smi_mols[0]))

    return products



if __name__ == "__main__":

    def args():

        import argparse
    
        parser = argparse.ArgumentParser()
        parser.add_argument("-reactants", help="", nargs="+", type=str, required=True)
        parser.add_argument("-reaction", help="", type=str, required=True)
        parser.add_argument("-out", help="", type=str, required=True)

        args = parser.parse_args()

        return args

    args = args()


    products = reactor(args.reactants, args.reaction)
    products = [" ".join(row) for row in products]

    with open(args.out, "w") as f:
        f.write("\n".join(products))
        f.close()

