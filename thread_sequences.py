#!/software/conda/envs/pyrosetta/

import os, sys
import pyrosetta
import numpy as np
import glob
import shutil
import argparse
import string
import random
import re
import subprocess
from pathlib import Path

##########################################################################################

from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.core.select.residue_selector import NeighborhoodResidueSelector
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta import ScoreFunction, create_score_function, get_fa_scorefxn
import pyrosetta.rosetta.core.scoring as scoring
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

##########################################################################################

parser = argparse.ArgumentParser(description='This script threads MPNN sequences onto a pdb backbone)')

parser.add_argument("--sym", type=str, required=True, help="symmetry")

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
sym = args.sym

##########################################################################################

#Parse MPNN .fa outputs and returns dict of name and seq
def parse_fasta(fasta, sym):
    with open(fasta, "r") as f:
        count = 0
        identifier = None
        sequence = []
        identifiers = []
        comp_seq = []
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                sequence = []
                count += 1
                line1 = line.split()
                score = [name for name in line1 if 'score' in name][0].rstrip(",")
                identifier = Path(fasta).stem + "_" + str(count) + "_" + score.replace("=", "_")
                sequence = []
                identifiers.append(identifier)
            else:
                sequence.append(line)
                seq = sequence[0].split("/")
                uniq_chains_ax_rot = set(seq)
                ABCDE_seq =  list(uniq_chains_ax_rot)[0]
                complete_sequence = int(sym) * ABCDE_seq
                comp_seq.append(complete_sequence)
        dict_id_seq = dict(zip(identifiers, comp_seq))
    return dict_id_seq

#Thread MPNN sequences onto a PDB
def thread_sequence(dict_seq, pdb_to_mutate):
    #tf = TaskFactory()
    #tf.push_back(operation.InitializeFromCommandline())
    #tf.push_back(operation.RestrictToRepacking())
    #packer = pack_min.PackRotamersMover()
    #packer.task_factory(tf)
    
    for name, seq in dict_seq.items():
        pose_ax_and_rot=pyrosetta.pose_from_pdb(pdb_to_mutate)
        mutated_pose = pyrosetta.rosetta.protocols.simple_moves.SimpleThreadingMover()
        mutated_pose.set_sequence(seq,1)
        mutated_pose.apply(pose_ax_and_rot)    
        #Pack the combined pose
        #packer.apply(pose_ax_and_rot)

        #Dumping axle and rotor alone
        print("Now dumping designed PDB... \n"  + str(pose_ax_and_rot) )
        pose_ax_and_rot.dump_pdb('mpnned_%s.pdb'%(name))

#############################################

pyrosetta.init('-nstruct 1 ') #-out:mute all  -ex1 -ex2 -beta

#############################################

def main():                      
    pdb_files = glob.glob( "sym*.pdb", recursive = True)
    print("About to thread all these pdbs: " + str(pdb_files))
    
    for pdb in pdb_files:
        #parse MPNN output sequences
        path_to_fasta = "output/temp_0.1/seqs/" + os.path.splitext(pdb)[0] + ".fa"
        parsed_seqs = parse_fasta(path_to_fasta, sym )
        
        #Thread new sequence onto PDB
        thread_sequence(parsed_seqs, pdb)
            
    print ("All done!")
    
main()
