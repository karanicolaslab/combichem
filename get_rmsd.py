import argparse
import glob, pandas as pd
from multiprocessing import Pool
from utils.io import read_sdf
from utils.rmsd import substructure_rmsd

def args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-r1", required=True)
    parser.add_argument("-r2", required=True)
    parser.add_argument("-core", required=True)
    parser.add_argument("-output", required=True)

    args = parser.parse_args()

    return args


def substructure__rmsd(ref, target, substructure):
    try:
        ref_mol = read_sdf(ref, add_hs=False)[0]
        target_mol = read_sdf(target, add_hs=False)[0]

        sub_mol = read_sdf(substructure, add_hs=False)[0]

        rmsd = substructure_rmsd(ref_mol, target_mol, sub_mol)

        return ref, target, rmsd[0], rmsd[1]
    except:
        return ref, target, None, None

if __name__ == "__main__":

    args = args()

    refs = glob.glob(args.r1 + "/*.sdf")
    targets = glob.glob(args.r2  + "/*.sdf")

    core = args.core

    data = []

    for r in refs:
        for t in targets:
            data.append([r,t,core])

    with Pool(64) as pool:
        data = pool.starmap(substructure__rmsd, data)

    data = pd.DataFrame(data, columns=["Ref", "Target", "RMSD", "Core"])
    data.to_csv(args.output, index=False)
