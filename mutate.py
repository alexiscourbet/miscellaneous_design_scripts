#!/usr/bin/python
import os
import re
import numpy as np
import pyrosetta
from math import *
import argparse

#############################################

###Give the script a pdb with chain A and B corresponding to each component.
###Give a fasta file with a bunch of sequences of A/B

#############################################

#init Rosetta
pyrosetta.init('-ex1 -ex2aro -nstruct 1 -beta') #-out:mute all

from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta import ScoreFunction, create_score_function, get_fa_scorefxn
import pyrosetta.rosetta.core.scoring as scoring

#############################################

argparser.add_argument("--pdb", type=str, help="pdb to mutate")
argparser.add_argument("--fasta", type=str, help="sequences to putate pdb into")
#argparser.add_argument("--suffix", type=str, help="suffix to add to name")
argparser.add_argument("--symA", type=str, help="symmetry to apply to pose A")
argparser.add_argument("--symB", type=str, help="symmetry to apply to pose B")

args = argparser.parse_args()
#############################################

def mutate(target_pose,sequence):
	for res in target_pose.

	pyrosetta.mutate_residue(target_pose, res, 'E')
              
	for i in range(len(sequence)):
	    pyrosetta.mutate_residue(target_pose, %d, 'E' %(i, str1[i]))
	yield(pose)

def read_fasta(fp):
        name, seqA, seqB = None, [], []
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name: 
                	yield (name, ''.join(seqA), ''.join(seqB))
                name, seqA, seqB = line, [], []
            else:
            	re.split(r'\t+', line)
                seqA.append(line[0])
                seqB.append(line[1])
        if name: 
            re.sub("\s+", "_", name.strip())
			re.sub(",", "_", name.strip())
        	yield (name, ''.join(seqA),''.join(seqB))

def repack(pose):
	scorefxn = pyrosetta.create_score_function("beta")
	task = standard_packer_task(pose) 
	task.restrict_to_repacking() 
	pack_mover = PackRotamersMover(scorefxn, task) 
	pack_mover.apply(pose) 
	yield(pose)

def symmetrize(pose,symmetry):
	SetupForSymmetryMover(symmetry).apply(pose)
	yield(pose)

###function to extract one chain from a pose
def pose_from_chain(pose, chain):
    r = [] # residues to keep
    for res in pose.residues:
        num = res.seqpos()
        chainletter = pose.pdb_info().chain(num)
        if chainletter == chain:
            r.append(num)
    return (Pose(pose, r[0], r[-1]))

#############################################

pose_to_mutate = pyrosetta.pose_from_pdb(args.pdb)
chainA_to_mutate = pose_from_chain(pose_to_mutate, 'A')
chainB_to_mutate = pose_from_chain(pose_to_mutate, 'B')

with open(args.fasta) as fp:
    for name, seqA, seqB in read_fasta(fp):
        print("Now working on..." + name, seqA, seqB)
		print ("Mutating pdbs according to fasta file")
		make_mutationA = mutate(chainA_to_mutate,seqA)
		make_mutationB = mutate(chainB_to_mutate,seqB)
		print ("Symmetrizing mutated pdbs...")
		symmetrize(make_mutationA,args.symA)
		symmetrize(make_mutationB,args.symB)
		combined_pose=make_mutationA.clone()
    	combined_pose.append_pose_by_jump(make_mutationB.clone(),1)
    	print ("Repaking...")
		repack(combined_pose)
		print ("Outputting mutated pdbs...")
		combined_pose.dump_pdb('mutated_%s.pdb'%(name))
		print ("Done")
