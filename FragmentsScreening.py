import os, sys, glob

if __name__ == "__main__":

	directory = sys.argv[1]
	sdfs = glob.glob(f"{directory}/*.sdf")
	done = [f"{directory}/" + f.split("/")[-1].split("__")[0].split("mini_")[-1] + ".sdf" for f in glob.glob(f"{directory}_minimization/*pdb")]

	sdfs = set(sdfs).difference(set(done))

	os.system(f"mkdir -p {directory}_mol2 {directory}_pdb_params {directory}_replaced_params {directory}_complexes {directory}_minimization")

	ref_pdb = "../mini_295_3BLQA_LG1_0001.pdb"
	ref_params = "../mini_295_3BLQA_LG1.params"

	ref_prot = "3BLQ_protein.pdb"
	klifs = "23 24 25 26 27 28 29 30 31 32 33 34 35 45 46 47 48 49 50 62 63 64 65 66 67 68 69 70 71 72 73 74 76 77 78 79 80 81 82 83 84 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 165 166 167 168 169 170 171"

	for sdf in sdfs:
		filename = sdf.split("/")[-1].split(".")[0]

		mol2 = f"{directory}_mol2/{filename}.mol2"
		pdb_params = f"{directory}_pdb_params/{filename}"

		am1bcc = f"assigncharges.py -method am1bcc -in {sdf} -out {mol2}"
		genpot = f"python ~/combichem/combichem/external_apps/generic_potencial/mol2genparams.py -s {mol2} --prefix={pdb_params}"

		replaced_params = f"{directory}_replaced_params/{filename}.params"
		replace = f"python ~/combichem/combichem/processing/Params.py -ref_pdb {ref_pdb} -ref_params {ref_params} -target_pdb {pdb_params}_0001.pdb -target_params {pdb_params}.params -out {replaced_params}"

		pdb_complex = f"{directory}_complexes/{filename}__{ref_prot}"
		concat = f"cat {ref_prot} {pdb_params}_0001.pdb > {pdb_complex}"

		min_pdb_complex = f"{directory}_minimization/mini_{filename}__{ref_prot}"
		min_log_complex = f"{directory}_minimization/mini_{filename}__{ref_prot}".replace(".pdb", ".log")

		minimization = f"python ~/combichem/combichem/processing/Screening.py -pdb {pdb_complex} -params {replaced_params} -output {min_pdb_complex} -residues {klifs} -beta True > {min_log_complex}"

		print("; ".join([am1bcc, genpot, replace, concat, minimization]))
