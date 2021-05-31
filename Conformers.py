"""Conformers generator based on RDKit engine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from Alignment import align
import utils.io as io

def generate_conformers(mol, num_confs, num_threads, substrucure=None,
                        max_iterations=1000, prune_rms_thresh=0.5):
    """Generate conformers for input molecule

    Args:
        mol (rdkit.Chem.rdchem.Mol): Input molecule
        num_confs (int): Max number of conformers for molecule
        num_threads (int): Number of threads for sampling
        substrucure (rdkit.Chem.rdchem.Mol, optional): Substructure to align
        max_iterations (int): TBD
        prune_rms_thresh (float): Diffence between conformers

    Returns:
        rdkit.Chem.rdchem.Mol: Molecule with conformers

    Raises:
        Exception: if mol is none, throw exception
    """
    params = AllChem.ETKDG()
    params.maxIterations = max_iterations
    params.numThreads = num_threads
    params.useBasicKnowledge = True
    params.useExpTorsionAnglePrefs = True
    params.pruneRmsThresh = prune_rms_thresh

    ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

    for i in ids:
        AllChem.MMFFOptimizeMolecule(mol, confId=i)

    Chem.rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol, numThreads=num_threads, maxIters=max_iterations, mmffVariant="MMFF94s")

    if substrucure is None:
        AllChem.AlignMolConformers(mol)
    else:
        mol = align(mol, substrucure)[0]

    if mol is None:
        raise Exception("Conformers cannot be generated")

    return mol



def main():
    import argparse
    """Summary

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-inp",
                        help="input smiles file",
                        type=str,
                        required=True)
    parser.add_argument("-out",
                        help="output SDF file",
                        type=str,
                        required=True)
    parser.add_argument("-num_conf",
                        help="Number of conformers for one conformers",
                        type=int,
                        default=100)
    parser.add_argument("-cpu",
                        help="Number of threads",
                        type=int,
                        default=1)
    parser.add_argument("-substructure",
                        help="Substucture to be fixed",
                        type=str,
                        default=None)
    parser.add_argument("-max_iterations",
                        help="",
                        type=int,
                        default=1000)
    parser.add_argument("-prune_rms_thresh",
                        help="",
                        type=float,
                        default=0.5)


    args = parser.parse_args()

    inp = args.inp
    out = args.out
    num_conf = args.num_conf
    cpu = args.cpu
    substructure = args.substructure
    max_iterations = args.max_iterations
    prune_rms_thresh = args.prune_rms_thresh

    # parse smiles
    with open(inp, "r") as smiles:
        smiles = [l for l in smiles.read().split("\n") if l != ""]

    smiles_mol = []

    for i,smi_row in enumerate(smiles):
        if " " in smi_row:
            smi, name = smi_row.split(" ")
        else:
            smi = smi_row
            name = "Mol_" + str(i+1)
        smi = io.read_smi(smi)
        smi.SetProp("_Name", name)
        smiles_mol.append(smi)

    if substructure is not None:
        substructure = io.read_sdf(substructure, add_hs=False)[0]

    for i, mol in enumerate(smiles_mol):
        mol = generate_conformers(mol, num_conf, cpu, substructure, max_iterations, prune_rms_thresh)

        if i == 0:
            io.write_sdf(mol, out, kekulize=False, status="w")
        else:
            io.write_sdf(mol, out, kekulize=False, status="a")


if __name__ == "__main__":
    main()
