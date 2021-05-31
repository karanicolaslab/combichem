from .mcs import substucture_search
from .io import get_atoms_names

def calculate_rmsd(mol1, mol2, atom_map=None):
    rmsd = 0

    for i, j in atom_map:
        mol1_coord = mol1.GetConformer(0).GetAtomPosition(i)
        mol2_coord = mol2.GetConformer(0).GetAtomPosition(j)
        dx = (mol2_coord.x - mol1_coord.x)**2
        dy = (mol2_coord.y - mol1_coord.y)**2
        dz = (mol2_coord.z - mol1_coord.z)**2
        rmsd += dx + dy + dz

    rmsd = rmsd / len(atom_map)
    rmsd = rmsd ** 0.5

    return rmsd


def find_closest_mcs(mol1, mol2, atom_maps):

    best_match = None
    best_rmsd = None

    for i, matches in enumerate(atom_maps):
        rmsd = calculate_rmsd(mol1, mol2, atom_map=matches)
        if best_rmsd is None or rmsd < best_rmsd:
            best_rmsd = rmsd
            best_match = i

    return atom_maps[best_match], best_rmsd


def substructure_rmsd(mol1, mol2, substructure):

    mol1_names = get_atoms_names(mol1)
    mol2_names = get_atoms_names(mol2)

    matches1 = [list(zip(*m)) for m in substucture_search(substructure, mol1)]
    matches2 = [list(zip(*m)) for m in substucture_search(substructure, mol2)]

    matches1 = {i:j for i,j in find_closest_mcs(substructure, mol1, matches1)[0]}
    matches2 = {i:j for i,j in find_closest_mcs(substructure, mol2, matches2)[0]}

    
    ref_keys = set(matches1.keys()).intersection(set(matches2.keys()))

    true_match = []
    for key in ref_keys:
        true_match.append((matches1[key], matches2[key]))

    best_match, best_rmsd  = find_closest_mcs(mol1, mol2, [true_match])

    # print({mol1_names[key]:mol2_names[value] for key, value in best_match})

    return best_rmsd, len(best_match)
