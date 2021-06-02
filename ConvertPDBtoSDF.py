import sys, os, glob
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

def args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-inp", required=True)
    parser.add_argument("-params", required=True)
    parser.add_argument("-out", required=True)
    
    args = parser.parse_args()

    return args

def reset_aromaticity(mol, idx):

    atom = mol.GetAtomWithIdx(idx)

    atom.SetIsAromatic(False)

    for bond in atom.GetBonds():
        bond.SetIsAromatic(False)
        if bond.GetBondType() is Chem.rdchem.BondType.AROMATIC:
            bond.SetBondType(Chem.rdchem.BondType.SINGLE)

    return mol

def correct_charge(mol, atom_id):
    DEFAULTS = {"O":2, "N":3}

    atom = mol.GetAtomWithIdx(atom_id)

    symbol = atom.GetSymbol()
    if symbol in DEFAULTS:
        bonds = atom.GetBonds()
        num_bonds = len(bonds)
        valences = [bond.GetBondTypeAsDouble() for bond in bonds]
        charge = sum(valences) - DEFAULTS[symbol]
        atom.SetFormalCharge(int(charge))

    return mol

def solve_problem(mol):

    problems = Chem.rdmolops.DetectChemistryProblems(mol)

    for problem in problems:
        try:
            problem_atoms = [problem.GetAtomIdx()]
        except:
            problem_atoms = problem.GetAtomIndices()

        problem_type = problem.GetType()

        if problem_type == "AtomValenceException":
            for idx in problem_atoms:
                mol = correct_charge(mol, idx)

        if problem_type == "AtomKekulizeException":
            # print("AtomKekulizeException", problem_atoms)
            for idx in problem_atoms:
                mol = reset_aromaticity(mol, idx)

        if problem_type == "KekulizeException":
            for idx, _ in enumerate(mol.GetAtoms()):
                mol = correct_charge(mol, idx)

            #print(write_smi(mol, kekuleSmiles=False))

            # print(problem_type, problem_atoms)
            # print(write_smi(mol))
            # print()

    if len(problems) != 0:
        mol = solve_problem(mol)

    Chem.SanitizeMol(mol)
    Chem.rdmolops.SetAromaticity(mol, model=Chem.rdmolops.AromaticityModel.AROMATICITY_MDL)
    Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)

    return mol



def parse_bonds(filename):

    bond_types = { "1":Chem.BondType.SINGLE,
                   "2":Chem.BondType.DOUBLE,
                   "3":Chem.BondType.TRIPLE,
                   "4":Chem.BondType.AROMATIC }

    bonds = []

    with open(filename, "r") as f:
        f = f.read().split("\n")

    for row in f:
        if row.startswith("BOND_TYPE"):
            row = row.replace(" RING", "")
            row = [elem for elem in row.split(" ")[1:] if elem != ""]
            row = [row[0], row[1], row[-1].replace("#ORGBND", "")]
            row[-1] = bond_types[row[-1]]
            bonds.append(row)

    return bonds

def read_pdb(filename, add_hs=True, remove_hs=False, proximity_bonding=False, sanitize=False, params=None):

    mol = open(filename, "r").read()
    mol = mol[0:mol.find("CONECT")]

    mol = Chem.MolFromPDBBlock(mol, removeHs=remove_hs,
                              proximityBonding=proximity_bonding,
                              sanitize=sanitize)

    if remove_hs:
        mol = remove_hydrogens(mol)

    atom_names = {}

    for i, atom in enumerate(mol.GetAtoms()):
        atom_name = atom.GetPDBResidueInfo().GetName()
        atom_names[atom_name.replace(" ","")] = i

    if params is not None:
        bonds = parse_bonds(params)

        mol_editable = Chem.RWMol(mol)

        for bond in bonds:
            if bond[0] in atom_names and bond[1] in atom_names:
                begin_atom_name = bond[0]
                end_atom_name = bond[1]
                bond_type = bond[2]

                begin_idx = atom_names[begin_atom_name]
                end_idx = atom_names[end_atom_name]

                mol_editable.AddBond(begin_idx, end_idx, order=bond_type)

        mol = mol_editable.GetMol()
        mol = solve_problem(mol)

        for i, _ in enumerate(mol.GetAtoms()):
            mol = correct_charge(mol, i)

        # print(write_smi(mol))
        Chem.SanitizeMol(mol)

    for i, atom in enumerate(mol.GetAtoms()):
        # print(i, mol.GetAtomWithIdx(i).GetPDBResidueInfo().GetName().strip())
        atom.UpdatePropertyCache()

    mol.SetProp("_Name", filename.split("/")[-1].split(".")[0])

    return mol


def write_sdf(mol, filename, kekulize=False, status="w"):

    with open(filename, status) as f:
        try:
            mol_name = mol.GetProp("_Name")
        except:
            mol_name = filename.split("/")[-1].split(".")[0]

        for i, conf in enumerate(mol.GetConformers()):
            mol.SetProp("_Name", mol_name + "_" + str(i))
            try:
                mb = Chem.MolToMolBlock(mol, kekulize=kekulize, confId=i)
            except:
                mb = Chem.MolToMolBlock(mol, kekulize=False, confId=i)
            f.write(mb.__str__() + "$$$$\n")
        f.close()


if __name__ == "__main__":

    args = args()

    mol = read_pdb(filename=args.inp, params=args.params, add_hs=False, proximity_bonding=False)
    write_sdf(mol, args.out, kekulize=True, status="w")




