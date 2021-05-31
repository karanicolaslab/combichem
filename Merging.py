from rdkit import Chem
from utils.rmsd import find_closest_mcs
from utils.mcs import substucture_search
import utils.io as io

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def delete_placeholder_hs(mol, core_atoms):
    need_delete = []

    for i, atom in enumerate(mol.GetAtoms()):
        if i not in core_atoms:
            neighboors = [neighboor.GetIdx() for neighboor in atom.GetNeighbors() if neighboor.GetIdx() not in core_atoms]
            if len(neighboors) == 0:
                if atom.GetSymbol() == "H":
                    need_delete.append(i)

    mol_editable = Chem.RWMol(mol)

    # print(need_delete, core_atoms)

    for atom_id in sorted(need_delete, reverse=True):
        mol_editable.RemoveAtom(atom_id)

    mol = mol_editable.GetMol()

    for delete_atom in sorted(need_delete, reverse=True):
        core_atoms = [elem - 1 if elem > delete_atom else elem for elem in core_atoms]

    return mol, core_atoms


def find_common_substructure(mol1, mol2, submol):

    matches_mol1 = [list(zip(*m)) for m in substucture_search(submol, mol1)]
    matches_mol2 = [list(zip(*m)) for m in substucture_search(submol, mol2)]

    matches_mol1_mol2 = []

    for m1 in matches_mol1:
        for m2 in matches_mol2:
            m1_dict = {i:j for i,j in m1}
            m2_dict = {i:j for i,j in m2}

            match_mol1_mol2 = set(m1_dict.keys()).intersection(set(m2_dict.keys()))
            match_mol1_mol2 = [(m1_dict[submol_atm_id], m2_dict[submol_atm_id]) for submol_atm_id in match_mol1_mol2]
            matches_mol1_mol2.append(match_mol1_mol2)

    matches_mol1_mol2, best_rmsd  = find_closest_mcs(mol1, mol2, matches_mol1_mol2)

    # matches_mol1 = {i:j for i,j in find_closest_mcs(submol, mol1, matches_mol1)[0]}
    # matches_mol2 = {i:j for i,j in find_closest_mcs(submol, mol2, matches_mol2)[0]}

    # matches_mol1_mol2 = set(matches_mol1.keys()).intersection(set(matches_mol2.keys()))

    # matches_mol1_mol2 = [(matches_mol1[submol_atm_id], matches_mol2[submol_atm_id]) for submol_atm_id in matches_mol1_mol2]
    # matches_mol1_mol2, best_rmsd  = find_closest_mcs(mol1, mol2, [matches_mol1_mol2])

    matches_mol1_mol2 = {i:j for i, j in matches_mol1_mol2}

    return matches_mol1_mol2


def merge(ref_pose, target_pose, substructure_pose):
    """Summary
    """
    ref_names = io.get_atoms_names(ref_pose)
    target_names = io.get_atoms_names(target_pose)

    matches = find_common_substructure(ref_pose, target_pose, substructure_pose)

    # print({ref_names[i]:target_names[j] for i,j in matches.items()})

    ref_pose, matches_ref = delete_placeholder_hs(ref_pose, list(matches.keys()))
    target_pose, matches_target = delete_placeholder_hs(target_pose, list(matches.values()))
    matches = dict(zip(matches_target, matches_ref))

    # print(io.write_smi(ref_pose))
    # print(io.write_smi(target_pose))

    ref_pose_editable = Chem.RWMol(ref_pose)

    added = {}

    for bond in target_pose.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()

        if begin_idx not in matches.keys() and begin_idx not in added.keys():
            begin_idx_new = ref_pose_editable.AddAtom(bond.GetBeginAtom())
            coord = target_pose.GetConformer(0).GetAtomPosition(begin_idx)
            ref_pose_editable.GetConformer(0).SetAtomPosition(begin_idx_new, coord)
            added[begin_idx] = begin_idx_new

        if end_idx not in matches.keys() and end_idx not in added.keys():
            end_idx_new = ref_pose_editable.AddAtom(bond.GetEndAtom())
            coord = target_pose.GetConformer(0).GetAtomPosition(end_idx)
            ref_pose_editable.GetConformer(0).SetAtomPosition(end_idx_new, coord)
            added[end_idx] = end_idx_new

        if begin_idx not in matches.keys() and end_idx not in matches.keys():
            ref_pose_editable.AddBond(added[begin_idx], added[end_idx], order=bond.GetBondType())

        if begin_idx in matches.keys() and end_idx not in matches.keys():
            ref_pose_editable.AddBond(matches[begin_idx], added[end_idx], order=bond.GetBondType())

        if begin_idx not in matches.keys() and end_idx in matches.keys():
            ref_pose_editable.AddBond(added[begin_idx], matches[end_idx], order=bond.GetBondType())

    ref_pose = ref_pose_editable.GetMol()

    # print(io.write_smi(ref_pose))

    ref_pose.UpdatePropertyCache()

    Chem.SanitizeMol(ref_pose)
    Chem.rdmolops.SetAromaticity(ref_pose, model=Chem.rdmolops.AromaticityModel.AROMATICITY_MDL)
    Chem.rdmolops.Kekulize(ref_pose, clearAromaticFlags=True)

    if ref_pose is None:
        raise Exception("Merging was unsuccessful!")

    return ref_pose




def main():
    import argparse
    """Summary
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-ref_pose", required=True)
    parser.add_argument("-target_pose", nargs="+", required=True)
    parser.add_argument("-substructure", required=True)
    parser.add_argument("-output", required=True)

    args = parser.parse_args()

    ref_pose = args.ref_pose
    target_pose = args.target_pose
    substructure = args.substructure
    output = args.output

    merged_mol = io.read_sdf(ref_pose, add_hs=False, remove_hs=False)[0]
    substructure = io.read_sdf(substructure, add_hs=False, remove_hs=False)[0]

    for pose in target_pose:
        pose = io.read_sdf(pose, add_hs=False, remove_hs=False)[0]
        merged_mol = merge(merged_mol, pose, substructure)

    io.write_sdf(merged_mol, output, kekulize=False, status="w")

if __name__ == "__main__":
    main()