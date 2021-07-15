import argparse
import sys, os, glob
from utils.io import read_sdf, write_sdf
from Alignment import align
from rdkit import Chem

def make_dirs(prefix, target_r_groups):

    if prefix != "":
        prefix += "_"

    os.system(f"mkdir -p {prefix}refs")
    os.system(f"mkdir -p {prefix}merged")
    os.system(f"mkdir -p {prefix}pdb_params")
    os.system(f"mkdir -p {prefix}replaced_params_ref")
    os.system(f"mkdir -p {prefix}complexes")
    os.system(f"mkdir -p {prefix}minimizations")

    for i in range(target_r_groups):
        os.system(f"mkdir -p {prefix}target_R{i+2}")
        os.system(f"mkdir -p {prefix}alns_R{i+2}")
        os.system(f"mkdir -p {prefix}replaced_params_R{i+2}")

    return f"{prefix}refs", [f"{prefix}target_R{i+2}" for i in range(target_r_groups)],\
           [f"{prefix}alns_R{i+2}" for i in range(target_r_groups)],\
           f"{prefix}merged", f"{prefix}pdb_params", f"{prefix}replaced_params_ref",\
           [f"{prefix}replaced_params_R{i+2}" for i in range(target_r_groups)],\
           f"{prefix}complexes", f"{prefix}minimizations"

def args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-out_prefix", default="")
    parser.add_argument("-ref_pdb", required=True)
    parser.add_argument("-ref_params", required=True)
    parser.add_argument("-target_pdb", nargs="+", required=True)
    parser.add_argument("-target_params", nargs="+", required=True)
    parser.add_argument("-substructure", required=True)
    parser.add_argument("-reaction", required=True)
    parser.add_argument("-chain", required=True)
    parser.add_argument("-resn", required=True)
    parser.add_argument("-klifs_seq", nargs="+", required=True)

    args = parser.parse_args()

    return args

def alignment(ref, targ, substructure, output):
    ref_pose = read_sdf(ref, add_hs=False)[0]
    target_pose = read_sdf(targ, add_hs=False)[0]

    substructure_pose = read_sdf(substructure, add_hs=False, remove_hs=True)[0]

    target_pose_aln, rmsd = align(target_pose, ref_pose, substructure=substructure_pose)
    write_sdf(target_pose_aln, output, status="w")

    return output, rmsd

def merging(ref, targets, core, output):
    p = os.path.realpath(__file__)
    p = p[0:p.rfind("/")+1]

    targets = " ".join(targets)

    cmd = f"python {p}/Merging.py -ref_pose {ref} -target_pose {targets} -substructure {core} -output {output}"
    os.system(cmd)

    mol = Chem.SDMolSupplier(output)[0]
    mol = Chem.RemoveHs(mol)
    # mol = Chem.AddHs(mol)

    with open(output, "w") as fwr:
        fwr.write(Chem.MolToMolBlock(mol))
        fwr.close()


def cmd_am1bcc_charge(input_sdf, output_mol):
    if os.getenv("OECHARGE"):
        os.system(f"python $OECHARGE -method am1bcc -in {input_sdf} -out {output_mol}")
    elif os.getenv("OBABEL"):
        os.system("obabel -isdf {input_sdf} -omol2 > {output_mol}")
    else:
        raise("Please install OEChem or OpenBabel")

def cmd_params(input_mol2, protein, reaction, chain, resn, pdb_params_prefix, complex_prefix):

    p = os.path.realpath(__file__)
    p = p[0:p.rfind("/")+1]

    os.system(f"python {p}/ParamsNCAACov.py -lig_mol {input_mol2} -protein {protein} \
                                        -out_pdb_params {pdb_params_prefix} \
                                        -reaction {reaction} -out_complex {complex_prefix} \
                                        -chain {chain} -resn {resn}")


def cmd_replace_params(ref_pdb, ref_params, input_pdb, input_params, output_params):
    p = os.path.realpath(__file__)
    p = p[0:p.rfind("/")+1]
    os.system(f"python {p}/Params.py -ref_pdb {ref_pdb} -ref_params {ref_params} -target_pdb {input_pdb} -target_params {input_params} -out {output_params}")

def minimization(ligand, protein, params, hinge):

    complexx = f"complexes/{ligand.split('/')[-1].split('.')[0]}.pdb"
    cmd = f"cat {protein} {ligand} > {complexx}"
    os.system(cmd)

    mini_complexx = f"minimizations/mini_{complexx.split('/')[-1]}"
    mini_log = f"minimizations/mini_{complexx.split('/')[-1].split('.')[0]}.log"
    cmd = f"python Screening.py -pdb {complexx} -params {params} -output {mini_complexx} -residues {hinge} > {mini_log}"

    os.system(cmd)


def extract_lig_prot(input_pdb, output_ligand_pdb, output_protein_pdb):
    with open(input_pdb, "r") as f:
        f = f.read().split("\n")

        atoms = []
        hetatms = []

        for line in f:
            if line == "":
                continue
            if line.startswith("ATOM"):
                atoms.append(line)
            elif line.startswith("HETATM") and "LG1" in line:
                hetatms.append(line)

        with open(output_ligand_pdb, "w") as f:
            f.write("\n".join(hetatms))
            f.close()

        with open(output_protein_pdb, "w") as f:
            f.write("\n".join(atoms))
            f.write("\nTER\n")
            f.close()


def cmd_convert_pdb_to_sdf(input_pdb, params, output_sdf):
    p = os.path.realpath(__file__)
    p = p[0:p.rfind("/")+1]

    os.system(f"python {p}/ConvertPDBtoSDF.py -inp {input_pdb} -params {params} -out {output_sdf}")


def cmd_concat(protein, input_pdb, output_pdb):
    os.system(f"cat {protein} {input_pdb} > {output_pdb}")


def cmd_minimization(input_pdb, input_params, input_klifs, output_pdb, output_log):
    p = os.path.realpath(__file__)
    p = p[0:p.rfind("/")+1]

    input_klifs = " ".join(input_klifs)

    os.system(f"python {p}/Screening.py -pdb {input_pdb} -params {input_params} -residues {input_klifs} -output {output_pdb} -no_bou_ubo > {output_log}")


if __name__ == "__main__":

    args = args()
    dirs = make_dirs(args.out_prefix, len(args.target_pdb))
    refs_dir = dirs[0]
    target_dirs = dirs[1]
    alns_dirs = dirs[2]
    merged_dir, pdb_params_dir, replaced_params_ref_dir = dirs[3:6]
    replaced_params_target_dirs = dirs[6]
    complexes_dir, minimization_dir = dirs[7::]

    pair_prefix = [args.ref_pdb] + args.target_pdb
    pair_prefix = [e.split("/")[-1].split(".")[0] for e in pair_prefix]
    pair_prefix = "_".join(pair_prefix)

    # extract ligand and proteins for reference
    input_pdb = args.ref_pdb
    ref_filename = input_pdb.split("/")[-1].split(".")[0].split("__")[0].split("mini_")[-1]
    output_ligand_pdb = f"{refs_dir}/{pair_prefix}_lig.pdb"
    output_protein_pdb = f"{refs_dir}/{pair_prefix}_prot.pdb"
    extract_lig_prot(input_pdb, output_ligand_pdb, output_protein_pdb)

    # convert ligands to sdf for reference
    input_pdb = output_ligand_pdb
    output_sdf = f"{refs_dir}/{pair_prefix}_lig.sdf"
    cmd_convert_pdb_to_sdf(input_pdb, args.ref_params, output_sdf)

    input_target_sdfs = []

    for i, tpdb in enumerate(args.target_pdb):
        # extract ligand and proteins for target
        input_pdb = tpdb
        filename = input_pdb.split("/")[-1].split(".")[0].split("__")[0].split("mini_")[-1]
        output_ligand_pdb = f"{target_dirs[i]}/{pair_prefix}_lig.pdb"
        output_protein_pdb = f"{target_dirs[i]}/{pair_prefix}_prot.pdb"
        extract_lig_prot(input_pdb, output_ligand_pdb, output_protein_pdb)

        # convert ligands to sdf for target
        input_pdb = output_ligand_pdb
        output_sdf = f"{target_dirs[i]}/{pair_prefix}_lig.sdf"
        cmd_convert_pdb_to_sdf(input_pdb, args.target_params[i], output_sdf)

        # align ligands to reference
        input_ref_sdf = f"{refs_dir}/{pair_prefix}_lig.sdf"
        input_target_sdf = f"{target_dirs[i]}/{pair_prefix}_lig.sdf"
        output_sdf = f"{alns_dirs[i]}/{pair_prefix}.sdf"

        alignment(input_ref_sdf, input_target_sdf, args.substructure, output_sdf)
        input_target_sdfs.append(input_target_sdf)

    # merging to one compound
    input_ref_sdf = input_ref_sdf
    input_ref_sdf_filename = input_ref_sdf.split("/")[-1].split("_lig")[0]
    input_target_sdfs = [f"alns_R{i+2}/{pair_prefix}.sdf" for i,f in enumerate(input_target_sdfs)]
    output_sdf = f"{merged_dir}/{pair_prefix}.sdf"
    merging(input_ref_sdf, input_target_sdfs, args.substructure, output_sdf)

    # charge assignment
    input_sdf = output_sdf
    compound_filename = output_sdf.split("/")[-1].split(".")[0]
    output_mol2 = f"{merged_dir}/{compound_filename}.mol2"
    cmd_am1bcc_charge(input_sdf, output_mol2)

    # Rosetta's param file generation
    input_mol2 = output_mol2
    output_prefix = f"{pdb_params_dir}/{compound_filename}"
    # cmd_params(input_mol2, output_prefix)
    cmd_params(input_mol2, f"{refs_dir}/{pair_prefix}_prot.pdb", args.reaction, args.chain, args.resn, f"{pdb_params_dir}/", f"{complexes_dir}/{compound_filename}.pdb")

    # Replace params for reference
    input_ref_pdb = f"{refs_dir}/{pair_prefix}_lig.pdb"
    input_ref_params = args.ref_params
    input_pdb = f"{output_prefix}_0001.pdb"
    input_params = f"{output_prefix}.params"
    output_params = f"{replaced_params_ref_dir}/{compound_filename}.params"
    cmd_replace_params(input_ref_pdb, input_ref_params, input_pdb, input_params, output_params)

    for i, tpdb in enumerate(args.target_pdb):
        # replace params for targets
        input_ref_pdb = f"{target_dirs[i]}/{pair_prefix}_lig.pdb"
        input_ref_params = args.target_params[i]

        input_pdb = f"{output_prefix}_0001.pdb"
        input_params = output_params

        output_params = f"{replaced_params_target_dirs[i]}/{compound_filename}.params"
        cmd_replace_params(input_ref_pdb, input_ref_params, input_pdb, input_params, output_params)

    # concat ligand and target protein structure
    output_pdb = f"{complexes_dir}/{compound_filename}.pdb"
    # cmd_concat(protein, input_pdb, output_pdb)

    # minimization of protein-ligand complex
    input_pdb = output_pdb
    input_params = output_params
    input_klifs = args.klifs_seq
    output_pdb = f"{minimization_dir}/mini_{compound_filename}.pdb"
    output_log = f"{minimization_dir}/mini_{compound_filename}.log"
    cmd_minimization(input_pdb, input_params, input_klifs, output_pdb, output_log)
