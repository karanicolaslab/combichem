import argparse, sys, os
from rdkit import Chem

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from MCS.mcs import substrucuture_matching

def AssignBondOrdersFromTemplate(refmol, mol):
    """ Updated version of original AssignBondOrdersFromTemplate

    file: <PYTHON_DIR>/site-packages/rdkit/Chem/AllChem.py

    """

    refmol2 = Chem.Mol(refmol)
    mol2 = Chem.Mol(mol)

    matches = substrucuture_matching(mol2, refmol2, remove_hs=False)[0]
    matching = {matches[1][i]:matches[0][i] for i in range(len(matches[0]))}

    # print(len(refmol2.GetAtoms()), len(mol2.GetAtoms()), len(matching))

    for b in refmol.GetBonds():
        atom1 = matching[b.GetBeginAtomIdx()]
        atom2 = matching[b.GetEndAtomIdx()]
        b2 = mol2.GetBondBetweenAtoms(atom1, atom2)
        b2.SetBondType(b.GetBondType())
        b2.SetIsAromatic(b.GetIsAromatic())


    for a in refmol.GetAtoms():
        a2 = mol2.GetAtomWithIdx(matching[a.GetIdx()])
        a2.SetHybridization(a.GetHybridization())
        a2.SetIsAromatic(a.GetIsAromatic())
        a2.SetNumExplicitHs(a.GetNumExplicitHs())
        a2.SetFormalCharge(a.GetFormalCharge())

    Chem.SanitizeMol(mol2)
    # Chem.SanitizeMol(mol2, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)

    return mol2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input_pdb", required=True)
    parser.add_argument("-input_sdf", required=True)
    parser.add_argument("-output", required=True)

    args = parser.parse_args()

    input_pdb = args.input_pdb
    input_sdf = args.input_sdf
    outpt = args.output

    mol_pdb = Chem.MolFromPDBFile(input_pdb, removeHs=False, proximityBonding=False)
    mol_sdf = Chem.SDMolSupplier(input_sdf, removeHs=False)[0]
    mol_sdf = Chem.AddHs(mol_sdf)

    if mol_pdb is None:
        print("PDB error: ", input_pdb)

    if mol_sdf is None:
        print("SDF error: ", input_pdb)


    print(mol_pdb.GetNumAtoms(), mol_sdf.GetNumAtoms(),)


    try:
        mol_pdb = AssignBondOrdersFromTemplate(mol_sdf, mol_pdb)
        mb = Chem.MolToMolBlock(mol_pdb, kekulize=False, includeStereo=True)

        if mol_pdb is None:
            print(input_pdb, input_sdf, "have some problems", sep=" ")
            return 0

        with open(outpt, "w") as f:
            f.write(mb.__str__())
            f.close()
    except Exception as e:
        print("Matching error: ", input_pdb)

if __name__ == "__main__":
    main()