import argparse
import os, sys, glob

def make_dirs(prefix):

    if prefix != "":
        prefix += "_"

    os.system(f"mkdir -p {prefix}sdf")
    os.system(f"mkdir -p {prefix}mol2")
    os.system(f"mkdir -p {prefix}pdb_params")
    os.system(f"mkdir -p {prefix}replaced_params")
    os.system(f"mkdir -p {prefix}complexes")
    os.system(f"mkdir -p {prefix}minimization")

    return f"{prefix}sdf", f"{prefix}mol2", f"{prefix}pdb_params", f"{prefix}replaced_params", \
           f"{prefix}complexes", f"{prefix}minimization"

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-out_prefix", default="")
    parser.add_argument("-smi", required=True)
    parser.add_argument("-num_confs", default=10)
    parser.add_argument("-substructure", required=True)
    parser.add_argument("-ref_pdb", required=True)
    parser.add_argument("-ref_param", required=True)
    parser.add_argument("-protein", required=True)
    parser.add_argument("-klifs_seq", nargs="+", required=True)
    
    args = parser.parse_args()

    return args


def cmd_conformers(input_smi, output_sdf, substructure, num_confs):

    if os.getenv("OEOMEGA"):
        # run omega
        os.system(f"$OEOMEGA -in {input_smi} -out {output_sdf} -commentEnergy true -prefix {output_sdf.split('.')[-1]} -sdEnergy true -warts true -maxconfs {num_confs} -strictstereo false")
    else:
        # run rdkit
        p = os.path.realpath(__file__)
        p = p[0:p.rfind("/")+1] 

        os.system(f"python {p}/Conformers.py -inp {input_smi} -out {output_sdf} -substructure {substructure} -num_conf {num_confs}")

def cmd_split(input_sdf, output_prefix):
    with open(input_sdf, "r") as f:
        f = f.read().split("$$$$\n")

        for i, mol in enumerate(f):
            if mol == "":
                continue

            name = mol.split("\n")[0]
            mol = mol + "$$$$\n"

            with open(f"{output_prefix}/{name}.sdf", "w") as fwr:
                fwr.write(mol)
                fwr.close()


def cmd_am1bcc_charge(input_sdf, output_mol):
    if os.getenv("OECHARGE"):
        os.system(f"python $OECHARGE -method am1bcc -in {input_sdf} -out {output_mol}")
    elif os.getenv("OBABEL"):
        os.system("obabel -isdf {input_sdf} -omol2 > {output_mol}")
    else:
        raise("Please install OEChem or OpenBabel")

def cmd_params(input_mol2, output_prefix):
    os.system(f"python $MOL2GENPARAMS -s {input_mol2} --prefix={output_prefix}")

def cmd_replace_params(ref_pdb, ref_params, input_pdb, input_params, output_params):
    p = os.path.realpath(__file__)
    p = p[0:p.rfind("/")+1] 
    os.system(f"python {p}/Params.py -ref_pdb {ref_pdb} -ref_params {ref_params} -target_pdb {input_pdb} -target_params {input_params} -out {output_params}")

def cmd_concat(protein, input_pdb, output_pdb):
    os.system(f"cat {protein} {input_pdb} > {output_pdb}")

def cmd_minimization(input_pdb, input_params, input_klifs, output_pdb, output_log):
    p = os.path.realpath(__file__)
    p = p[0:p.rfind("/")+1] 

    input_klifs = " ".join(input_klifs)

    os.system(f"python {p}/Screening.py -pdb {input_pdb} -params {input_params} -residues {input_klifs}  -output {output_pdb} > {output_log}")

if __name__ == "__main__":

    args = args()
    sdf_dir, mol2_dir, pdb_params_dir, replaced_params_dir, complexes_dir, minimization_dir = make_dirs(args.out_prefix)

    filename = args.smi.split("/")[-1].split(".")[0]

    # generate conformers
    input_smi = args.smi
    num_confs = args.num_confs
    substructure = args.substructure
    if args.out_prefix != "":
        output_sdf = f"{args.out_prefix}/{filename}.sdf"
    else:
        output_sdf = f"{filename}.sdf"
    cmd_conformers(input_smi, output_sdf, substructure, num_confs)

    # split conformers file on single structures
    input_sdf = output_sdf
    output_prefix = f"{sdf_dir}/"
    cmd_split(input_sdf, output_prefix)

    # process each conformer
    for sdf in glob.glob(f"{output_prefix}/*sdf"):
        filename = sdf.split("/")[-1].split(".")[0]

        # charge assignment
        input_sdf = sdf
        output_mol2 = f"{mol2_dir}/{filename}.mol2"
        cmd_am1bcc_charge(input_sdf, output_mol2)

        # Rosetta's param file generation
        input_mol2 = output_mol2
        output_prefix = f"{pdb_params_dir}/{filename}"
        cmd_params(input_mol2, output_prefix)

        # replace params 
        ref_pdb = args.ref_pdb
        ref_params = args.ref_param
        input_pdb = f"{output_prefix}_0001.pdb"
        input_params = f"{output_prefix}.params"
        output_params = f"{replaced_params_dir}/{filename}.params"
        cmd_replace_params(ref_pdb, ref_params, input_pdb, input_params, output_params)

        # concat ligand and target protein structure
        protein = args.protein
        protein_name = protein.split("/")[-1].split(".")[0]
        input_pdb = input_pdb
        output_pdb = f"{complexes_dir}/{filename}__{protein_name}.pdb"
        cmd_concat(protein, input_pdb, output_pdb)

        # minimization of protein-ligand complex
        input_pdb = output_pdb
        input_params = output_params
        input_klifs = args.klifs_seq
        output_pdb = f"{minimization_dir}/mini_{filename}__{protein_name}.pdb"
        output_log = f"{minimization_dir}/mini_{filename}__{protein_name}.log"
        cmd_minimization(input_pdb, input_params, input_klifs, output_pdb, output_log)
        

