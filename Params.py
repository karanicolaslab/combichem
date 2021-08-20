import utils.io as io
from utils.mcs import substucture_search
from utils.rmsd import find_closest_mcs


def get_charges(file_r):
    """Extract partial charges from params file

    Args:
        file (str): Filename of params file

    Returns:
        dict: Dictionary of partial charges
    """
    charges_dict = {}

    for line in file_r:
        if line[0:4] == "ATOM":
            line = line.split(" ")
            line = [elem for elem in line if elem != ""]
            charges_dict[line[1]] = [line[2], line[-1]]

    return charges_dict


def replace_charge(reference, ref_params, target, target_params):
    """Summary

    Args:
        target (str): Filename of target params file
        source (str): Filename of target params file
        matches_dict (dict): dictionary of matches of target and source

    Returns:
        str: Params files with updated options
    """
    # print(ref_pdb, target_pdb, sep="\t", end="\t")

    reference_names = io.get_atoms_names(reference)
    target_names = io.get_atoms_names(target)

    matched_atoms = [list(zip(*m)) for m in substucture_search(target, reference)]

    matched_atoms, rmsd = find_closest_mcs(target, reference, matched_atoms)
    matched_atoms = {target_names[i]:reference_names[j] for i,j in matched_atoms}

    # print(matched_atoms, len(matched_atoms), sep="\t", end="\n")
    # print(matched_atoms)

    source_charges = get_charges(ref_params)

    pool = []

    for line in target_params:
        if line[0:4] == "ATOM":
            line = line.split(" ")
            line = [elem for elem in line if elem != ""]
            if line[1] in matched_atoms:
                line[2] = source_charges[matched_atoms[line[1]]][0]
                line[4] = source_charges[matched_atoms[line[1]]][1]
                line = "ATOM %-4s %-4s %-4s %6.3f" % (
                    line[1], line[2], line[3], float(line[4]))
            else:
                line = "ATOM %-4s %-4s %-4s %6.3f" % (
                    line[1], line[2], line[3], float(line[4]))
            pool.append(line)
        elif line[0:10] == "NBR_RADIUS":
            pool.append("NBR_RADIUS 999.0000")
        else:
            pool.append(line)

    return "\n".join(pool), len(matched_atoms)



def save_file(params, out):
    with open(out, "w") as f:
        f.write(params)
        f.close()

def repl_charge(ref_pdb, ref_params, target_pdb, target_params, out):
    reference = io.read_pdb(ref_pdb, params=ref_params)
    target = io.read_pdb(target_pdb, params=target_params)
    
    replaced_params, num = replace_charge(reference, target)
    
    if replaced_params is not None:
        save_file(replaced_params, out)

    return num



def main():
    import argparse
    """Summary
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-ref_pdb",
                        help="Reference pdb",
                        type=str,
                        required=True)
    parser.add_argument("-ref_params",
                        help="Reference params file",
                        type=str,
                        required=True)
    parser.add_argument("-target_pdb",
                        help="Target pdb",
                        type=str,
                        required=True)
    parser.add_argument("-target_params",
                        help="Target params file",
                        type=str,
                        required=True)
    parser.add_argument("-out",
                        help="Output folder",
                        type=str,
                        required=True)

    args = parser.parse_args()

    ref_pdb = args.ref_pdb
    ref_params = args.ref_params
    target_pdb = args.target_pdb
    target_params = args.target_params
    out = args.out

    repl_charge(ref_pdb, ref_params, target_pdb, target_params, out)


if __name__ == "__main__":
    main()


