from rdkit import Chem
from rdkit.Chem import rdMolAlign

from utils.mcs import substucture_search
from utils.io import read_sdf, write_sdf
from utils.rmsd import find_closest_mcs

def align(target_pose, ref_pose, substructure=None):
    """Aligmnent of target_pose against reference

    Args:
        ref_pose (TYPE): Description
        target_pose (TYPE): Description

    Returns:
        TYPE: Description
    """
    if substructure is None:
        substructure = Chem.Mol(ref_pose)

    matches1 = [list(zip(*m)) for m in substucture_search(substructure, ref_pose)]
    matches2 = [list(zip(*m)) for m in substucture_search(substructure, target_pose)]

    matches1 = {i:j for i,j in find_closest_mcs(substructure, ref_pose, matches1)[0]}
    matches2 = {i:j for i,j in find_closest_mcs(substructure, target_pose, matches2)[0]}

    ref_keys = set(matches1.keys()).intersection(set(matches2.keys()))

    true_match = []
    for key in ref_keys:
        true_match.append((matches2[key], matches1[key]))

    best_match, best_rmsd  = find_closest_mcs(target_pose, ref_pose, [true_match])

    rmsd = rdMolAlign.AlignMol(target_pose, ref_pose, atomMap=best_match)
    rdMolAlign.AlignMolConformers(target_pose, atomIds=[i for i, k in true_match])

    return target_pose, rmsd




def main():
    import argparse
    """ Align target molecule against reference
    
    Usage: alignment.py -ref <sdf> -target <sdf> -substructure <sdf> -out <sdf>
    
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-ref", help="", type=str, required=True)
    parser.add_argument("-target", help="", type=str, required=True)
    parser.add_argument("-substructure", help="", type=str, default=None)
    parser.add_argument("-out", help="", type=str, required=True)

    args = parser.parse_args()

    ref_pose_name = args.ref
    target_pose_name = args.target

    if args.substructure is not None:
        substructure = args.substructure
    else:
        substructure = ref_pose_name

    out_name = args.out

    ref_pose = read_sdf(ref_pose_name, add_hs=False)[0]
    target_pose = read_sdf(target_pose_name, add_hs=False)

    substructure_pose = read_sdf(substructure, add_hs=False)[0]

    for i, pose in enumerate(target_pose):
        target_pose_aligned = align(pose, ref_pose, substructure_pose)[0]
        write_sdf(target_pose_aligned, out_name, status="w")

if __name__ == "__main__":
    main()
