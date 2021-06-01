import os, sys, glob

def make_dirs(prefix=""):

    if prefix != "":
        prefix += "_"

    os.system(f"mkdir -p {prefix}mol2")
    os.system(f"mkdir -p {prefix}pdb_params")
    os.system(f"mkdir -p {prefix}replaced_params")
    os.system(f"mkdir -p {prefix}complexes")
    os.system(f"mkdir -p {prefix}minimization")

    return f"{prefix}mol2", f"{prefix}pdb_params", f"{prefix}replaced_params",
           f"{prefix}complexes", f"{prefix}minimization"

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-sdf", required=True)
    parser.add_argument("-prot", required=True)
    parser.add_argument("-ref_charge_pdb", required=True)
    parser.add_argument("-ref_charge_param", required=True)
    parser.add_argument("-klifs-seq", nargs="+", required=True)
    parser.add_argument("-out_prefix", default="")
    
    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = args()
    mol2_dir, pdb_params_dir, replaced_params_dir, complexes_dir, minimization_dir = make_dirs(args["out_prefix"])

    filename = args["sdf"].split("/")[-1].split(".")[0]

    mol2 = f"{mol2_dir}/{filename}.mol2"
    pdb_params = f"{pdb_params_dir}/{filename}"

    am1bcc = f"assigncharges.py -method am1bcc -in {sdf} -out {mol2}"
    genpot = f"python external_apps/generic_potencial/mol2genparams.py -s {mol2} --prefix={pdb_params}"

    replaced_params = f"{replaced_params_dir}/{filename}.params"
    replace = f"python Params.py -ref_pdb {ref_pdb} -ref_params {ref_params} -target_pdb {pdb_params}_0001.pdb -target_params {pdb_params}.params -out {replaced_params}"

    pdb_complex = f"{complexes_dir}/{filename}__{ref_prot}"
    concat = f"cat {ref_prot} {pdb_params}_0001.pdb > {pdb_complex}"

    min_pdb_complex = f"{minimization_dir}/mini_{filename}__{ref_prot}"
    min_log_complex = f"{minimization_dir}/mini_{filename}__{ref_prot}".replace(".pdb", ".log")

    minimization = f"python Screening.py -pdb {pdb_complex} -params {replaced_params} -output {min_pdb_complex} -residues {klifs} -beta True > {min_log_complex}"

    print("; ".join([am1bcc, genpot, replace, concat, minimization]))
