import argparse, os
from itertools import permutations
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_to_params import Params

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-lig_mol2", help="", type=str, required=True)
    parser.add_argument("-reaction", help="", type=str, required=True)
    parser.add_argument("-protein", help="", type=str, required=True)
    parser.add_argument("-out_pdb_params", help="", type=str, required=True)
    parser.add_argument("-out_complex", help="", type=str, required=True)
    parser.add_argument("-chain", help="", type=str, default="A")
    parser.add_argument("-resn", help="", type=int, default=12)

    args = parser.parse_args()

    return args

def parse_params(file):

    def parse_atoms(line):
        name = [e for e in line.split(" ") if e != ""][1]
        return {name:line}


    atoms = {}

    file = open(file, "r").read().split("\n")
    for row in file:
        if row.startswith("ATOM"):
            atoms.update(parse_atoms(row))
            
    return atoms

def conjugate_mols(ligand, reaction):

    reaction = open(reaction,"r").read().split("\n")[0]
    rxn = AllChem.ReactionFromSmarts(reaction)

    reactants = [ligand]
    
    product = rxn.RunReactants(reactants)
    Chem.SanitizeMol(product[0][0])

    return product[0][0]

def extract_atoms_names(conjugate, ligand):

    ligand_names = {i:a.GetPropsAsDict()["_TriposAtomName"] for i,a in enumerate(ligand.GetAtoms())}

    conjugate_names = {} 

    for i, c_atom in enumerate(conjugate.GetConformer(0).GetPositions()):
        x1, y1, z1 = c_atom

        for j, l_atom in enumerate(ligand.GetConformer(0).GetPositions()):
            x2, y2, z2 = l_atom
            rmsd =  (x1 - x2) ** 2
            rmsd += (y1 - y2) ** 2
            rmsd += (z1 - z2) ** 2
            rmsd = rmsd ** 0.5

            if rmsd == 0.00:
                conjugate_names[i] = ligand_names[j]

    return conjugate_names


def replace_atom_names(params, labels):
    for i, a in enumerate(params.mol.GetAtoms()):
        atom_name = params._get_PDBInfo_atomname(a, throw=True)
        if atom_name == "LOWER" or atom_name == "UPPER":
            continue
        params.rename_atom(atom_name, f"Z{i}")    

    for i, a in enumerate(params.mol.GetAtoms()):
        atom_name = params._get_PDBInfo_atomname(a, throw=True)
        if atom_name == "LOWER" or atom_name == "UPPER":
            continue
        params.rename_atom(atom_name, labels[i])

    return params


def clean_params(params, atoms):
    
    rows = []
    for r in params.split("\n"):
        if r.startswith("ATOM"):
            name = [e for e in r.split(" ") if e != ""][1]
            r = atoms[name]
        if r.startswith("BOND"):
            if r[:5] == "BOND ":
                r = r.replace("BOND", "BOND_TYPE") + " 1"
            r += " #ORGBND" + [e for e in r.split(" ") if e != ""][-1]

        rows.append(r)

    return "\n".join(rows)

def gen_params(mol2, output):
    os.system(f"python $MOL2GENPARAMS -s {mol2} --prefix={output}")

def concatenate_with_protein(ligand_pdb, protein_pdb, resn, output):
    with open(ligand_pdb, "r") as lig_atoms:
        lig_atoms = [r for r in lig_atoms.read().split("\n") if r.startswith("ATOM")]
    
    new_pdb = []

    with open(protein_pdb, "r") as protein_atoms:
        protein_atoms = protein_atoms.read().split("\n")

        flag_done = False

        for r in protein_atoms:
            if r.startswith("ATOM"):
                if not flag_done:
                    if int(r[22:26].strip()) - 1 == resn:
                        for rr in lig_atoms:
                            new_pdb.append(rr)
                        flag_done = True
            new_pdb.append(r)

    new_pdb = "\n".join(new_pdb)

    with open(output, "w") as fwr:
        fwr.write(new_pdb)
        fwr.close()

if __name__ == "__main__":
    
    args = args()

    output = args.out_pdb_params + "/" + args.lig_mol2.split("/")[-1].split(".")[0]
    gen_params(args.lig_mol2, output)

    lig_atoms = parse_params(output + ".params")

    ligand = Chem.MolFromMol2File(args.lig_mol2, sanitize=True, removeHs=False)

    t = Chem.MolFromSmarts("[H:6][#7-:5][C:3]([H:4])([#6:1]([H])=[O:2])[C:7]([H:8])([H:9])[#16:10]")

    conjugate = conjugate_mols(ligand, args.reaction)
    conjugate_names = extract_atoms_names(conjugate, ligand)

    conjugate_params = Params.from_mol(conjugate, name="LG1")
    conjugate_params.polish_mol(resi=args.resn, chain=args.chain)
    conj_params = conjugate_params.dumps().replace("\n", " \n")
    conj_pdb = Chem.MolToPDBBlock(conjugate_params.dummyless)
    conj_pdb = [e for e in conj_pdb.split("\n") if e.startswith("ATOM") or e.startswith("HETATM")]
    conj_pdb_left = [s.split("LG1")[0] for s in conj_pdb]
    conj_pdb_left = "\n".join(conj_pdb_left)
    conj_pdb_right = [s.split("LG1")[1] for s in conj_pdb]

    temp_ids = {}

    for i, atom in enumerate(conjugate_params.mol.GetAtoms()):
        oldname = atom.GetPDBResidueInfo().GetName()
        if "LOWER" in oldname or "UPPER" in oldname:
            continue

        newname = f"Z{i+1}"
        conj_params = conj_params.replace(f" {oldname.strip()} ", f" {newname} ")
        conj_pdb_left = conj_pdb_left.replace(f" {oldname} ", f"  {newname} " + (3 - len(newname.strip())) * " ")
        temp_ids[i] = newname


    for i, atom in enumerate(conjugate_params.mol.GetAtoms()):
        if i not in temp_ids:
            continue
        
        oldname = temp_ids[i]
        newname = conjugate_names[i]
        newname = newname + " " * (len(oldname) - len(newname))

        conj_params = conj_params.replace(f" {oldname} ", f" {newname} ")
        conj_pdb_left = conj_pdb_left.replace(f" {oldname} ", f" {newname} ")

    conj_pdb = []
    for i, lr in enumerate(conj_pdb_left.split("\n")):
        row = lr + "LG1" + conj_pdb_right[i]
        row = [e for e in row.split(" ") if e != ""]

        atom = row[0]
        seq_atom = int(row[1])
        atomn = row[2]
        resn = row[3]
        chain = row[4]
        resi = int(row[5])
        x = float(row[6])
        y = float(row[7])
        z = float(row[8])
        temp = float(row[9])
        occ = float(row[10])
        atomtype = row[11]

        _a = "{:6s}{:5d} {:<4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}"
        _a += "{:6.2f}{:6.2f}     {:<4s}{:<2s}{:2s}"
        row = _a.format(atom, seq_atom, atomn, "", resn, chain, resi, "", x, y, z, temp, occ, "", "", atomtype)
        conj_pdb.append(row)

    conj_pdb = "\n".join(conj_pdb)

    os.system(f"rm -rf {output}_0001.pdb {output}.params")    

    with open(f"{output}.params", "w") as fwr:
        conj_params = clean_params(conj_params, lig_atoms)
        fwr.write(conj_params)
        fwr.close()

    with open(f"{output}_0001.pdb", "w") as fwr:
        fwr.write(conj_pdb)
        fwr.close()

    concatenate_with_protein(f"{output}_0001.pdb", args.protein, args.resn, args.out_complex)
