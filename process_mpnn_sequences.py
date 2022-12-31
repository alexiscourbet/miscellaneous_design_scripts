#!/software/conda/envs/pyrosetta/

import os, sys
import numpy as np
import glob
import shutil
import argparse
import string
import random
import re
from pathlib import Path

##########################################################################################


##########################################################################################

parser = argparse.ArgumentParser(description='This script parses MPNN outputs sequences and concatenates them input AF2 parsable fasta')   

parser.add_argument("--number_of_chains", type=str, required=True, help="number of chains ")

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
number_of_chains = args.number_of_chains

##########################################################################################

#Parse MPNN .fa outputs and returns dict of name and seq
def parse_fasta(fasta, number_of_chains):
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
                A_seq = seq[0]
                complete_sequence = ( int(number_of_chains) - 1 ) * ( A_seq + 'UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU' ) +  A_seq
                comp_seq.append(complete_sequence)
        dict_id_seq = dict(zip(identifiers, comp_seq))
    return dict_id_seq

#############################################

#############################################

def main():                      
    pdb_files = glob.glob( "*.pdb", recursive = True)
    print("About to process all these pdbs sequences: " + str(pdb_files))
    full_dict = {}

    for pdb in pdb_files:
        print("Now working on: " + str(pdb))

        #parse MPNN output sequences
        path_to_fasta = "output/temp_0.1/seqs/" + os.path.splitext(pdb)[0] + ".fa"
        parsed_seqs = parse_fasta(path_to_fasta, number_of_chains )
        print(parsed_seqs)
        full_dict.update(parsed_seqs)

    f = open("sequences_to_fold.fasta", "w")
    for k in full_dict.keys():
        f.write(">" + str(k) + "\n" + str(full_dict[k]) + "\n")
    f.close()

    print ("All done!")
    
main()


		
