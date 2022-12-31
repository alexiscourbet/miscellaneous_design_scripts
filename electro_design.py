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

from rosetta.core.kinematics import MoveMap
from rosetta.core.kinematics import FoldTree
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from rosetta.core.simple_metrics import metrics
from rosetta.core.select import residue_selector as selections
from rosetta.core import select
from rosetta.core.select.movemap import *

from rosetta.protocols import minimization_packing as pack_min
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
    

##########################################################################################

parser = argparse.ArgumentParser(description='This script finds interface between machine components (axle and rotor) and designs them using only charged electrostatic interactiosn accorss the interface. Give the script pdbs corresponding to axle/rotor assembly (axle should start with chain A, and axle+rotor res numbering should be continuous)')

#parser.add_argument('--axle_and_rotor', default=None,required=False,help='axle component plus rotor component')
parser.add_argument("--sym_axle", type=int, required=True, help="symmetry of axle, if C5 -> 5, if D4 -> 4 ")
parser.add_argument("--sym_rotor", type=int, required=True, help="symmetry of rotor, if C4 -> 4, if D3 -> 3")
parser.add_argument('--nstructs', type=int, default=None, required=True, help='Number of designs to output')

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
sym_axle = args.sym_axle
sym_rotor = args.sym_rotor
nstructs = args.nstructs

##########################################################################################

#find the residue identities of the interface between axle and rotor
def find_the_interface_res(pdbs_axle_and_rotor,sym_ax,sym_rot):
    ABCDetc = re.sub(r'([A-Z])(?!$)', r'\1,', string.ascii_uppercase[0:int(sym_ax)])
    FGHIJKetc = re.sub(r'([A-Z])(?!$)', r'\1,', string.ascii_uppercase[int(sym_ax):int(sym_ax)+int(sym_rot)])
    chains_axle = ChainSelector(ABCDetc)
    chains_rotor = ChainSelector(FGHIJKetc)
    interface_on_ax = NeighborhoodResidueSelector(chains_rotor, 10.0, False)
    interface_on_rot = NeighborhoodResidueSelector(chains_axle, 10.0, False)
    interface_by_vector = pyrosetta.rosetta.core.select.residue_selector.InterGroupInterfaceByVectorSelector(interface_on_ax, interface_on_rot)
    interface_by_vector.cb_dist_cut(8.5) #8
    interface_by_vector.vector_angle_cut(75)
    interface_by_vector.vector_dist_cut(9.5) #9
    
    pose_axrot=pyrosetta.pose_from_pdb(pdbs_axle_and_rotor)
    by_vector = interface_by_vector.apply(pose_axrot)

    number_of_resi_total = pose_axrot.total_residue()
    #TO DES
    seq_list_to_des = []
    for sequence_position in range(1, number_of_resi_total):
        if by_vector[sequence_position] == True:
            seq_list_to_des.append(sequence_position)
        else:
            pass
    
    res_chain_end_id = list(pyrosetta.rosetta.core.pose.chain_end_res(pose_axrot))
    
    def count(list1, l, r):
        c = 0
        sublist = []
        # traverse in the list1
        for x in list1:
            # condition check
            if x - 1 >= l and x<= r:
                c+= 1
                sublist.append(x)
        return sublist
         
    number_of_chains = sym_ax + sym_rot
    
    renumbered_list_ax = []
    renumbered_list_ax.append(count(seq_list_to_des, 0, res_chain_end_id[0]))
    renumbered_list_rot = []
    renumbered_list_rot.append(count(seq_list_to_des, res_chain_end_id[sym_ax - 1], res_chain_end_id[sym_ax]))
    new_renumbered_list_ax = []
    list_temp_ax = []
    new_renumbered_list_rot = []
    list_temp_rot = []
    
    for i in range(0, sym_ax):
        length_chain_ax = res_chain_end_id[1] - res_chain_end_id[0]
        new_renumbered_list_ax = [[element + i * length_chain_ax for element in sublist] for sublist in renumbered_list_ax]
        list_temp_ax.append(new_renumbered_list_ax)
    for i in range(0, sym_rot):
        length_chain_rot = res_chain_end_id[-1] - res_chain_end_id[-2] 
        new_renumbered_list_rot = [[element + i * length_chain_rot for element in sublist] for sublist in renumbered_list_rot]
        list_temp_rot.append(new_renumbered_list_rot)
    
    flat_list_rot = [x for xs in list_temp_rot for x in xs]
    flat_list_rot2 =  [x for xs in flat_list_rot for x in xs]
    flat_list_ax = [x for xs in list_temp_ax for x in xs]
    flat_list_ax2 = [x for xs in flat_list_ax for x in xs]
    joined_lists = list_temp_ax + list_temp_rot
    import itertools
    joined_lists_n = list(itertools.chain.from_iterable(joined_lists))
    flat_list_n = [x for xs in joined_lists_n for x in xs]
    
    #NOT TO DES
    seq_list_not_to_des = []
    for sequence_position in range(1, number_of_resi_total):
        if by_vector[sequence_position] == False:
            seq_list_not_to_des.append(sequence_position)
        else:
            pass
    
    renumbered_list_ax_nodes = []
    renumbered_list_ax_nodes.append(count(seq_list_not_to_des, 0, res_chain_end_id[0]))
    renumbered_list_rot_nodes = []
    renumbered_list_rot_nodes.append(count(seq_list_not_to_des, res_chain_end_id[sym_ax - 1], res_chain_end_id[sym_ax]))
    new_renumbered_list_ax_nodes = []
    list_temp_ax_nodes = []
    new_renumbered_list_rot_nodes = []
    list_temp_rot_nodes = []
    
    for i in range(0, sym_ax):
        length_chain_ax_nodes = res_chain_end_id[1] - res_chain_end_id[0]
        new_renumbered_list_ax_nodes = [[element + i * length_chain_ax_nodes for element in sublist] for sublist in renumbered_list_ax_nodes]
        list_temp_ax_nodes.append(new_renumbered_list_ax_nodes)
    for i in range(0, sym_rot):
        length_chain_rot_nodes = res_chain_end_id[-1] - res_chain_end_id[-2] 
        new_renumbered_list_rot_nodes = [[element + i * length_chain_rot_nodes for element in sublist] for sublist in renumbered_list_rot_nodes]
        list_temp_rot_nodes.append(new_renumbered_list_rot_nodes)
    
    flat_list_rot_nodes = [x for xs in list_temp_rot for x in xs]
    flat_list_rot2_nodes =  [x for xs in flat_list_rot for x in xs]
    flat_list_ax_nodes = [x for xs in list_temp_ax for x in xs]
    flat_list_ax2_nodes = [x for xs in flat_list_ax for x in xs]
    joined_lists_nodes = list_temp_ax_nodes + list_temp_rot_nodes
    import itertools
    joined_lists_n_nodes = list(itertools.chain.from_iterable(joined_lists_nodes))
    flat_list_n_nodes = [x for xs in joined_lists_n_nodes for x in xs]
    print(flat_list_n_nodes)

    return flat_list_ax2, flat_list_rot2, flat_list_n_nodes

#run FD
def run_rosetta_design(list_ax, list_rot, list_no_design, sym_ax, sym_rot, pdb_):
    pose_ax_and_rot=pyrosetta.pose_from_pdb(pdb_)
    #Select residues on axle and rotor
    selector_ax = ResidueIndexSelector()
    for residue in list_ax:
        selector_ax.append_index(residue)       
    selector_rot = ResidueIndexSelector()
    for residue in list_rot:
        selector_rot.append_index(residue)
                
    #select residues to lock
    selector_axrot = ResidueIndexSelector()
    for residue in list_no_design:
        selector_axrot.append_index(residue)
                     
    # The task factory accepts all the task operations
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    
    # Disable design (repack only) on residues to lock that are not at the interface
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), selector_axrot))
        
    # Enable design on axle
    aa_axle_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    aa_axle_to_design.aas_to_keep("DEQTNSYA")
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        aa_axle_to_design, selector_ax))
        
    # Enable design on rotor
    aa_rotor_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    aa_rotor_to_design.aas_to_keep("KRQTNSYA")
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
            aa_rotor_to_design, selector_rot))
        
    # Convert the task factory into a PackerTask
    packer = pack_min.PackRotamersMover()
    packer.task_factory(tf)
    packer.apply(pose_ax_and_rot)
    pdb_name = os.path.splitext(pdb_)[0]
    pose_ax_and_rot.dump_pdb('Designed_%s.pdb'%(pdb_name))

#create folder for MPNN runs
def make_new_dir(pdb_name):
    dirName = "output_FD_" + os.path.splitext(pdb_name)[0]
    if not os.path.exists(dirName):
        os.mkdir(dirName)
        print("Directory " , dirName ,  " created ")
    else:    
        print("Directory " , dirName ,  " already exists")
        
#############################################

pyrosetta.init('-nstruct ' + str(nstructs) + ' -beta') #-out:mute all  -ex1 -ex2

#############################################

def main():                      
    pdb_files = glob.glob( "*.pdb", recursive = True)
    print("About to design all these pdbs: " + str(pdb_files))
    #create output folder
    
    for pdb in pdb_files:
        #make_new_dir(pdb)
        print("Now working on: " + str(pdb))
        #extract interface
        flat_ax2, flat_rot2, flat_nodes = find_the_interface_res(pdb,sym_axle,sym_rotor)
        #run FD
        run_rosetta_design(flat_ax2, flat_rot2, flat_nodes, sym_axle, sym_rotor, pdb)
        
    print ("All done!")
    
main()


		
