class PDB:
    def __init__(self, filename):
        self.atoms = self._read_file(filename)

    def _read_file(self, filename):
        with open(filename, "r") as pdb:
            pdb = pdb.read().split("\n")
            pdb = [row for row in pdb if row.startswith("HETATM")]

            return self._parse_lines(pdb)
            
    def _parse_lines(self, lines):
        atoms = {}

        for line in lines:
            name = line[12:16].strip()
            x = line[30:38].strip()
            y = line[38:46].strip()
            z = line[46:54].strip()

            atoms[name] = (float(x),float(y),float(z))

        return atoms


def fmcs(ref_pdb, target_pdb):

    matches = {}

    for ref_atom_name, ref_atom_coords in ref_pdb.atoms.items():

        lowest_rmsd = 999
        lowest_atom_name = "None"

        ref_atom_elem = "".join([char for char in ref_atom_name if char.isalpha()])
        for target_atom_name, target_atom_coords in target_pdb.atoms.items():
            target_atom_elem = "".join([char for char in target_atom_name if char.isalpha()])
            if ref_atom_elem == target_atom_elem:
                dx = (ref_atom_coords[0] - target_atom_coords[0]) ** 2
                dy = (ref_atom_coords[1] - target_atom_coords[1]) ** 2
                dz = (ref_atom_coords[2] - target_atom_coords[2]) ** 2
                rmsd = (dx + dy + dz) ** 0.5

                if rmsd < lowest_rmsd:
                    lowest_rmsd = rmsd
                    lowest_atom_name = target_atom_name

        if lowest_rmsd < 1.0:
            matches[ref_atom_name] = lowest_atom_name

    return matches


def get_charges(file):
    """Extract partial charges from params file
    Args:
        file (str): Filename of params file
    Returns:
        dict: Dictionary of partial charges
    """
    with open(file, "r") as file_r:
        file_r = file_r.read().split("\n")

    charges_dict = {}

    for line in file_r:
        if line[0:4] == "ATOM":
            line = line.split(" ")
            line = [elem for elem in line if elem != ""]
            charges_dict[line[1]] = [line[2], line[-1]]

    return charges_dict

def rmsd(target_pdb, ref_pdb, matches):
    rmsd = 0

    for target_name, ref_name in matches.items():
        target_atom_coords = target_pdb.atoms[target_name]
        ref_atom_coords = ref_pdb.atoms[ref_name]

        dx = (target_atom_coords[0] - ref_atom_coords[0]) ** 2
        dy = (target_atom_coords[1] - ref_atom_coords[1]) ** 2
        dz = (target_atom_coords[2] - ref_atom_coords[2]) ** 2

        rmsd += dx + dy + dz

    rmsd = (rmsd / len(matches)) ** 0.5
    return rmsd


def replace_charge(ref_pdb, ref_params, target_pdb, target_params):
    """Summary
    Args:
        target (str): Filename of target params file
        source (str): Filename of target params file
        matches_dict (dict): dictionary of matches of target and source
    Returns:
        str: Params files with updated options
    """
    print(ref_pdb, target_pdb, sep="\t", end="\t")

    ref_pdb = PDB(ref_pdb)
    target_pdb = PDB(target_pdb)

    matched_atoms = fmcs(target_pdb, ref_pdb)

    print(rmsd(target_pdb, ref_pdb, matched_atoms), end="\t")

    print(matched_atoms, len(matched_atoms), sep="\t", end="\n")

    source_charges = get_charges(ref_params)

    with open(target_params, "r") as file_r:
        file_r = file_r.read().split("\n")

    pool = []

    for line in file_r:
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

    return "\n".join(pool)



def save_file(params, out):
    with open(out, "w") as f:
        f.write(params)
        f.close()

def repl_charge(ref_pdb, ref_params, target_pdb, target_params, out):
    replaced_params = replace_charge(ref_pdb, ref_params, target_pdb, target_params)
    
    if replaced_params is not None:
        save_file(replaced_params, out)



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