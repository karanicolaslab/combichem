"""Input/output function for processing SDF, SMILE and PDB structures
"""
from rdkit import Chem

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


def cat_complex(pdb1, pdb2, output):
    
    with open(pdb1, "r") as pdb1_file:
        pdb1_file = pdb1_file.read().split("\n")

    with open(pdb2, "r") as pdb2_file:
        pdb2_file = pdb2_file.read().split("\n")

    with open(output, "w") as f_wr:
        f_wr.write("\n".join(pdb1_file))
        f_wr.write("\n".join(pdb2_file))
        f_wr.close()


def reset_elem(input, output):

    lines = []

    with open(input,"r") as file:
        file = file.read().split("\n")

    for line in file:
        atomName = line[12:15].strip()
        atomElem = line[76:78].strip()

        if atomElem not in atomName:
            elem = "".join([symbol for symbol in atomName if not symbol.isdigit()])
            elem = "{0:>2}".format(elem)
            line = list(line)
            line[76:78] = elem[0], elem[1]
            line = "".join(line)

        lines.append(line)


    with open(output, "w") as file:
        file.write("\n".join(lines))
        file.close()    

def divide_complex(pdb, ligand_filename, protein_filename):
    ligand_rows = []
    protein_rows = []

    with open(pdb, "r") as pdb_file:
        pdb_file = pdb_file.read().split("\n")

    for line in pdb_file:
        if line.startswith("ATOM"):
            protein_rows.append(line)
        elif line.startswith("HETATM"):
            ligand_rows.append(line)

    with open(ligand_filename, "w") as f:
        f.write("\n".join(ligand_rows))
        f.close()

    reset_elem(ligand_filename, ligand_filename)


    with open(protein_filename, "w") as f:
        f.write("\n".join(protein_rows))
        f.close()


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


def reset_aromaticity(mol, idx):

    atom = mol.GetAtomWithIdx(idx)

    atom.SetIsAromatic(False)
    
    for bond in atom.GetBonds():
        bond.SetIsAromatic(False)
        if bond.GetBondType() is Chem.rdchem.BondType.AROMATIC:
            bond.SetBondType(Chem.rdchem.BondType.SINGLE)

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

            print(write_smi(mol, kekuleSmiles=False))

            # print(problem_type, problem_atoms)
            # print(write_smi(mol))
            # print()

    if len(problems) != 0:
        mol = solve_problem(mol)

    Chem.SanitizeMol(mol)
    Chem.rdmolops.SetAromaticity(mol, model=Chem.rdmolops.AromaticityModel.AROMATICITY_MDL)
    Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
    
    return mol



def read_pdb(filename, add_hs=True, remove_hs=False, proximity_bonding=False, sanitize=False, params=None):

    with open(filename, "r") as pdb:
        pdb = pdb.read().split("\n")
        pdb = [p for p in pdb if p.startswith("ATOM") or p.startswith("HETATM")]
        pdb = "\n".join(pdb)
    
    mol = Chem.MolFromPDBBlock(filename, removeHs=remove_hs,
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

def remove_hydrogens(mol):
    """Summary

    Remove hydrogens from molecule

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit container of input molecule

    Returns:
        rdkit.Chem.rdchem.Mol: input molecule without hydrogens
    """
    mol = Chem.RemoveHs(mol)

    ids = [i for i, atom in enumerate(mol.GetAtoms()) if atom.GetSymbol() == "H"]

    mol_editable = Chem.RWMol(mol)
    
    for atom_id in sorted(ids, reverse=True):
        mol_editable.RemoveAtom(atom_id)
        
    mol = mol_editable.GetMol()

    return mol


def add_hydrogens(mol, add_coords=True):
    """Summary

    Add hydrogents into molecule. Warning: add_coords is not guarantee
    correct position of hydrogens

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit container of input molecule
        add_coords (bool, optional): Set coordinates of hydrogens based on input

    Returns:
        rdkit.Chem.rdchem.Mol: input molecule with hydrogens
    """
    return Chem.AddHs(mol, addCoords=add_coords)


def read_sdf(filename, add_hs=True, remove_hs=False):
    """Convert SDF into RDKit mol

    Args:
        filename (str): Name of File
        add_hs (bool, optional): Description
        remove_hs (bool, optional): Remove hydrogens in input molecule

    Returns:
        list: list of RDKit molecules

    Raises:
        Exception: If input mol is None, throw exception
    """
    mols = Chem.SDMolSupplier(filename, removeHs=remove_hs, sanitize=False)

    if mols is None:
        raise Exception("Cannot read {0} structure".format(filename))

    mols = [mol for mol in mols]

    if remove_hs:
        mols = [remove_hydrogens(mol) for mol in mols]    

    if add_hs:
        for i, mol in enumerate(mols):
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^
                             Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
            mols[i] = add_hydrogens(mol, add_coords=True)

    return mols


def read_smi(smi, add_hs=True, remove_hs=False):
    """Convert smiles into RDKit mol

    Args:
        smi (str): smiles string
        add_hs (bool, optional): Add hydrogens to input molecule
        remove_hs (bool, optional): Description

    Returns:
        rdkit.Chem.rdchem.Mol: RDKit container of input molecule

    Raises:
        Exception: If input mol is None, throw exception
    """
    parser_params = Chem.SmilesParserParams()
    parser_params.removeHs = remove_hs
    parser_params.sanitize = False

    mol = Chem.MolFromSmiles(smi, params=parser_params)
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^
                     Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)

    if mol is None:
        raise Exception("Cannot read {0}".format(smi))

    if add_hs:
        mol = add_hydrogens(mol, add_coords=True)

    return mol


def write_smi(mol, kekuleSmiles=True, remove_hs=False):
    """Convert smiles into RDKit mol

    Args:
        smi (str): smiles string
        add_hs (bool, optional): Add hydrogens to input molecule
        remove_hs (bool, optional): Description

    Returns:
        rdkit.Chem.rdchem.Mol: RDKit container of input molecule

    Raises:
        Exception: If input mol is None, throw exception
    """
    if remove_hs:
        mol = remove_hydrogens(mol)

    smi = Chem.MolToSmiles(mol, kekuleSmiles=kekuleSmiles)

    return smi



def write_sdf(mol, filename, kekulize=False, status="w"):
    """Write SDF into file

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit container molecule
        filename (str): Name of output file
        kekulize (bool, optional): Aromatic or Kekule form in rings
        status (str, optional): Overwrite or append data in file
    """
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


def write_sdf_as_text(structure, prefix="./"):
    name = prefix + "/" + structure.split("\n")[0] + ".sdf"

    with open(name, "w") as f:
        f.write(structure)
        f.close()


def split_sdf(filename):
    with open(filename, "r") as f:
        f = f.read().split("$$$$\n")

    f = [elem + "$$$$\n" for elem in f if elem != ""]

    return f


def get_atoms_names(mol):
    """Extract atom names from PDB

    Args:
        mol (rdkit.mol): Input molecule

    Returns:
        dict: Dictionary of pairs rdkit_id:pdb_atom_name
    """
    dic = {}

    for idx in range(mol.GetNumAtoms()):
        try:
            dic[idx] = mol.GetAtomWithIdx(idx).GetPDBResidueInfo().GetName().strip()
        except:
            pass

    return dic
