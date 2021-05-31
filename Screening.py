import pyrosetta
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.scoring import ScoreType, CA_rmsd
from pyrosetta.rosetta.core.scoring import get_score_function, all_atom_rmsd
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.protocols.simple_moves import SuperimposeMover
from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover
from pyrosetta.rosetta.protocols.grafting import delete_region
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.protocols.toolbox.pose_manipulation import superimpose_pose_on_subset_CA
from pyrosetta.rosetta.utility import vector1_unsigned_long
from pyrosetta.rosetta.std import list_unsigned_long_t

def init(params=None, beta=True, mute=["protocols.simple_moves.SuperimposeMover",
                                       "protocols.moves.RigidBodyMover",
                                       "core.pack.dunbrack.RotamerLibrary",
                                       "basic.io.database",
                                       "core.scoring.etable",
                                       "core.chemical.GlobalResidueTypeSet",
                                       "core.scoring.elec.util",
                                       "core.scoring.P_AA",
                                       "core.scoring.GenericBondedPotential",
                                       "core.scoring.GenericBondedEnergy"]):
    
    init_str = "-flip_HNQ -no_optH true -ignore_zero_occupancy false"

    if params is not None:
        init_str += " -extra_res_fa {0}".format(params)

    if mute is not None:
        init_str += " -mute {0}".format(" ".join(mute))

    if beta:
        init_str += " -beta"


    pyrosetta.init(init_str)


def get_residue_num(pose, resname):
    for i in range(1, pose.size() + 1):
        if pose.residue(i).name3() == resname:
            return i

def superimpose(target, reference):
    sp_mover = SuperimposeMover()
    sp_mover.set_reference_pose(reference, 1, reference.size() - 1)
    sp_mover.set_target_range(1, (target.size() - 1))
    sp_mover.apply(target)


def get_score_terms(pose, res_num, scorefxn):
    scorefxn(pose)
    
    score_type = ScoreType.fa_intra_rep

    scores = pose.energies().residue_total_energies(res_num).weighted_string_of( scorefxn.weights() )

    return scores

def get_total(pose, scorefxn):
    scorefxn(pose)
    score_type = ScoreType.total_score
    return pose.energies().total_energies()[score_type]

def get_interaction_energy(pose, scorefxn):
    ubo_ps = pose.clone()
    rb_jump = ubo_ps.num_jump()
    ubo_dist = 1000000.0
    trans_mover = RigidBodyTransMover(ubo_ps, rb_jump)
    trans_mover.step_size(ubo_dist)
    trans_mover.apply(ubo_ps)

    difference = get_total(pose, scorefxn) - get_total(ubo_ps, scorefxn)

    return difference


def minimize_complex(pose,
                     scorefxn,
                     min_func="lbfgs_armijo_nonmonotone",
                     tolerance=0.000001,
                     max_iter=2000):


    mm = MoveMap()
    mm.set_chi(True)
    mm.set_bb(True)
    mm.set_jump(True)

    mover = MinMover(movemap_in=mm,
                     scorefxn_in=scorefxn,
                     min_type_in=min_func,
                     tolerance_in=tolerance,
                     use_nb_list_in=True,
                     deriv_check_in=False,
                     deriv_check_verbose_in=False)

    mover.max_iter(max_iter)
    mover.apply(pose)


def save_pose(pose, scorefxn, filename):
    pose.dump_scored_pdb(filename, scorefxn)



def get_hbonds(pose):

    hbonds_return = []

    hbonds_list = pose.get_hbonds()

    for hbond in hbonds_list.hbonds():
        acc_atom_num = hbond.acc_atm()
        acc_res_num = hbond.acc_res()

        acc_atom_name = pose.residue(acc_res_num).atom_name(acc_atom_num)
        acc_res_name = pose.residue(acc_res_num).name3()
        acc_pdb_num = pose.pdb_info().number(acc_res_num)

        don_atom_num = hbond.don_hatm()
        don_res_num = hbond.don_res()

        don_atom_name = pose.residue(don_res_num).atom_name(don_atom_num)
        don_res_name = pose.residue(don_res_num).name3()
        don_pdb_num = pose.pdb_info().number(don_res_num)
    
        hbonds_return.append({"ACC_RES":acc_res_name,
                              "ACC_NUM":acc_pdb_num,
                              "ACC_ATOM":acc_atom_name.replace(" ",""),
                              "DON_RES":don_res_name,
                              "DON_NUM":don_pdb_num,
                              "DON_ATOM":don_atom_name.replace(" ","")})

    return hbonds_return

def get_coordinates(pose, res_id, atom_names):
    lst = []

    for atom in atom_names:
        lst.append(list(pose.residue(res_id).atom(atom).xyz()))

    return lst

def rmsd_score(coord1, coord2):
    rmsd = 0
    for i in range(len(coord1)):
        dx = (coord2[i][0] - coord1[i][0]) ** 2
        dy = (coord2[i][1] - coord1[i][1]) ** 2
        dz = (coord2[i][2] - coord1[i][2]) ** 2
        rmsd += dx + dy + dz
    rmsd = (rmsd / len(coord1)) ** 0.5

    return rmsd


def screening(pdb, output, lig_name="LG1", tolerance=0.000001, max_iter=2000, residues=None, substructure=None):
    """Summary

    Args:
        pdb (str): Filename of input PDB file
        params (str): Params file of ligand in PDB file
        output (str, optional): Output name of minimized pdb
        tolerance (float, optional): Convergence threshold of minimization
        max_iter (int, optional): Number of iterations during minimization
    """

    bound_pose = pose_from_pdb(pdb)
    non_minimized_pose = bound_pose.clone()

    scorefxn = get_score_function()

    minimize_complex(bound_pose, scorefxn, tolerance=tolerance)

    vector = vector1_unsigned_long()
    long_list = list_unsigned_long_t()

    rosetta_res_id = {bound_pose.pdb_info().number(i):i for i in range(1, bound_pose.size() + 1)}

    for i, res_num in enumerate(residues):
        vector.append(rosetta_res_id[res_num])
        long_list.append(rosetta_res_id[res_num])

    ca_rmsd = superimpose_pose_on_subset_CA(bound_pose, non_minimized_pose, vector)

    save_pose(bound_pose, scorefxn, output)

    # printing energy score
    print("Residues: " + " ".join([str(elem) for elem in residues]))

    lig_res_num = get_residue_num(non_minimized_pose, lig_name)
    protein_energy = get_total(non_minimized_pose, scorefxn)
    lig_score_terms = get_score_terms(non_minimized_pose, lig_res_num, scorefxn)
    

    print("Initial score of complex: ", protein_energy)
    print("Initial score terms of ligand:", lig_score_terms)

    lig_score_terms = get_score_terms(bound_pose, lig_res_num, scorefxn)
    protein_energy = get_total(bound_pose, scorefxn)

    print("Final score of complex: ", protein_energy)
    print("Final fa_intra_rep of ligand:", lig_score_terms)

    energy = get_interaction_energy(bound_pose, scorefxn)
    print("Energy of bound-unbound complexes: ", energy)

    # printing rmsd
    all_rmsd = all_atom_rmsd(bound_pose, non_minimized_pose, long_list)

    print("Non-minimized/Minimized ca-RMSD: ", ca_rmsd)
    print("Non-minimized/Minimized all-atom-RMSD:", all_rmsd)

    hbonds = get_hbonds(bound_pose)

    for i, bond in enumerate(hbonds):
        if bond["ACC_RES"] == lig_name or bond["DON_RES"] == lig_name:
            print("Acceptor-Donor HBonds: ", bond["ACC_RES"], bond["ACC_NUM"], bond["ACC_ATOM"], bond["DON_RES"], bond["DON_NUM"], bond["DON_ATOM"])

    ligand_res_id = list_unsigned_long_t()
    ligand_res_id.append(lig_res_num)

    ligand_rmsd = all_atom_rmsd(bound_pose, non_minimized_pose, ligand_res_id)
    print("RMSD of ligand substructure: ", ligand_rmsd)


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-pdb", help="", type=str, required=True)
    parser.add_argument("-params", help="", type=str, required=True)
    parser.add_argument("-output", help="", type=str, required=True)
    parser.add_argument("-residues", help="", nargs="+", type=int, default=None)
    parser.add_argument("-beta", help="", type=bool, default=True)

    args = parser.parse_args()

    pdb = args.pdb
    params = args.params
    output = args.output
    residues = args.residues

    init(params=params)
    screening(pdb, output, residues=residues)


if __name__ == "__main__":
    main()