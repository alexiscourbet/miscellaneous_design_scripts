#!/software/conda/envs/pyrosetta/

import argparse
import sys
import pyrosetta
import pyrosetta.distributed.tasks.rosetta_scripts
import pyrosetta.distributed.io
import pyrosetta.distributed.packed_pose
#Core Includes
from rosetta.core.kinematics import MoveMap
from rosetta.core.kinematics import FoldTree
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from rosetta.core.simple_metrics import metrics
from rosetta.core.select import residue_selector as selections
from rosetta.core import select
from rosetta.core.select.movemap import *
#Protocol Includes
from rosetta.protocols import minimization_packing as pack_min

####################################################################################################################

parser = argparse.ArgumentParser(description='This script thread sequences onto axle and rotor parts of rotary machine..')

parser.add_argument('--axle_pdb', 
					type=str,
                    default=None,
                    required=True,
                    help='Path to axle pdb')

parser.add_argument('--rotor_pdb',
                    type=str,
                    required=True,
                    default=None, 
                    help='Path to rotor pdb')  

parser.add_argument('--axle_seq', 
					type=str,
                    default=None,
                    required=True,
                    help='The new sequence to mutate the axle with')

parser.add_argument('--rotor_seq',
                    type=str,
                    required=True,
                    default=None, 
                    help='The new sequence to mutate the axle with')            

parser.add_argument('--axle_sym', 
					type=str,
                    default=None,
                    required=False,
                    help='Symmetry to apply to the axle. Example: D5')

parser.add_argument('--rotor_sym',
                    type=str,
                    required=False,
                    default=None, 
                    help='Symmetry to apply to the rotor. Example: C3')                                                                                                                                                                                    

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

axle = args.axle_pdb
rotor = args.rotor_pdb
new_axle_seq = args.axle_seq
new_rotor_seq = args.rotor_seq
sym_axle = args.axle_sym
sym_rotor = args.rotor_sym

####################################################################################################################
#init Rosetta
pyrosetta.init("-mute all")

tf = TaskFactory()
tf.push_back(operation.InitializeFromCommandline())
tf.push_back(operation.RestrictToRepacking())
packer = pack_min.PackRotamersMover()
packer.task_factory(tf)

#Thread new sequence
def thread_sequence(pose,seq):
    mutated_pose = pyrosetta.rosetta.protocols.simple_moves.SimpleThreadingMover()
    mutated_pose.set_sequence(seq,1)
    mutated_pose.apply(pose)

def make_sym(pose,sym):
    setup_for_symmetry = pyrosetta.rosetta.protocols.symmetry.SetupForSymmetryMover()
    setup_for_symmetry.process_symmdef_file('/home/acourbet/symm/' + sym + "_Z.sym")
    setup_for_symmetry.apply(pose)
    reset_sym = pyrosetta.rosetta.protocols.symmetry.ExtractAsymmetricPoseMover()
    reset_sym.apply(pose)

####################################################################################################################

pose_axle=pyrosetta.pose_from_pdb(axle)
pose_rotor=pyrosetta.pose_from_pdb(rotor)

axle_sequence = pose_axle.sequence()
rotor_sequence = pose_rotor.sequence()

print('mutating axle from ' + axle_sequence + ' to ' + new_axle_seq + ' with symmetry ' + sym_axle )
print('mutating rotor from ' + rotor_sequence + ' to ' + new_rotor_seq + ' with symmetry ' + sym_rotor )

#thread sequences onto ASUs
thread_sequence(pose_axle,new_axle_seq)
thread_sequence(pose_rotor,new_rotor_seq)

#generate symmetry of axle and rotor
make_sym(pose_axle,sym_axle)
make_sym(pose_rotor,sym_rotor)

#Dumping axle and rotor alone
print("Now dumping axle and rotor components alone... \n" + 'Axle = \n'  + str(pose_axle) + '\n' + 'Rotor = \n' + str(pose_rotor))
pose_axle.dump_pdb( axle + sym_axle + '_mutated.pdb')
pose_rotor.dump_pdb( rotor + sym_axle + '_mutated.pdb')

#combine poses
combined_clone=pose_axle.clone()
combined_clone.append_pose_by_jump(pose_rotor.clone(),1)

#Pack the combined pose
packer.apply(combined_clone)

#Dump the PDB
print("Now dumping axle and rotor components alone... \n" + 'Axle + Rotor  = \n'  + str(combined_clone))
combined_clone.dump_pdb( axle + rotor + sym_axle + sym_rotor + '_mutated.pdb')
