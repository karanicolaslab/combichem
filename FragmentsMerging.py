import sys, os, glob
import pandas as pd
from multiprocessing import Pool
from combichem.utils.io import read_sdf, write_sdf
from combichem.processing.Alignment import align

def alignment(ref, targ, output):
    ref_pose = read_sdf(ref, add_hs=False)[0]
    target_pose = read_sdf(targ, add_hs=False)[0]

    target_pose_aln, rmsd = align(target_pose, ref_pose)
    write_sdf(target_pose_aln, output, status="w")

    return output, rmsd

def merging(r1, r2, core, output):

    cmd = f"python /fccc/users/karanicolaslab/andriag/combichem/combichem/processing/Merging.py -ref_pose {r1} -target_pose {r2} -substructure {core} -output {output}"
    os.system(cmd)

def generate_params(sdf, output):

    cmd = f"convert.py {sdf} {output}.mol2"
    os.system(cmd)
    cmd = f"python /fccc/users/karanicolaslab/andriag/combichem/combichem/external_apps/generic_potencial/mol2genparams.py -s {output}.mol2 --prefix={output}"
    os.system(cmd)

def replace_charge(r_pdb, r_params, t_pdb, t_params, output):

    cmd = f"python /fccc/users/karanicolaslab/andriag/combichem/combichem/processing/Params.py -ref_pdb {r_pdb} -ref_params {r_params} -target_pdb {t_pdb} -target_params {t_params} -out {output}"
    os.system(cmd)


def minimization(ligand, protein, params, hinge):

    complexx = f"complexes/{ligand.split('/')[-1].split('.')[0]}.pdb"
    cmd = f"cat {protein} {ligand} > {complexx}"
    os.system(cmd)

    mini_complexx = f"minimizations/mini_{complexx.split('/')[-1]}"
    mini_log = f"minimizations/mini_{complexx.split('/')[-1].split('.')[0]}.log"
    cmd = f"python /fccc/users/karanicolaslab/andriag/combichem/combichem/processing/Screening.py -pdb {complexx} -params {params} -output {mini_complexx} -residues {hinge} > {mini_log}"

    os.system(cmd)


if __name__ == "__main__":

    os.system("mkdir -p R2_aln complexes Merged minimizations ReplacedParams_R1 ReplacedParams_R2")

    r1_file = "R1/mini_R1_C3v0_3_0__3BLQ_protein.sdf"
    r2_file = "R2/mini_R2_C1v0_1__3BLQ_protein.sdf"
    core = sys.argv[1]

    hinge = "23 24 25 26 27 28 29 30 31 32 33 34 35 45 46 47 48 49 50 62 63 64 65 66 67 68 69 70 71 72 73 74 76 77 78 79 80 81 82 83 84 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 165 166 167 168 169 170 171"

    r1 = r1_file.split("/mini_")[-1].split("__")[0]
    r2 = r2_file.split("/mini_")[-1].split("__")[0]
    
    r2_aln = f"R2_aln/{r1}_by_{r2}.sdf"
    #alignment(r1_file, r2_file, r2_aln)

    merged_compound = f"Merged/{r1}_{r2}.sdf"
    merging(r1_file, r2_aln, core, merged_compound)

    params_prefix = merged_compound.split(".")[0]
    generate_params(merged_compound, params_prefix)

    t_pdb = f"{params_prefix}_0001.pdb"
    t_params = f"{params_prefix}.params"

    r1_pdb = glob.glob(f"R1/{r1}_0001.pdb")[0]
    r1_params = glob.glob(f"R1/{r1}.params")[0]

    t_r1_params = f"ReplacedParams_R1/{t_params.split('/')[-1]}"
    replace_charge(r1_pdb, r1_params, t_pdb, t_params, t_r1_params)

    t_params = t_r1_params

    r2_pdb = glob.glob(f"R2/{r2}_0001.pdb")[0]
    r2_params = glob.glob(f"R2/{r2}.params")[0]

    t_r2_params = f"ReplacedParams_R2/{t_params.split('/')[-1]}"
    replace_charge(r2_pdb, r2_params, t_pdb, t_params, t_r2_params)

    protein = r1_file.replace("sdf","pdb")
    minimization(t_pdb, protein, t_r2_params, hinge)
