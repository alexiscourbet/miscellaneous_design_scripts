import pymol
from pymol import cmd
import glob, re
import os

##script to aling many PDBs onto target and output RMSD

file_list = glob.glob( "allpdbs*pdb")
target = "target.pdb"
extension = re.compile( '(^.*[\/]|\.(pdb|ent|brk))' )
target_name = extension.sub('',target)
file_list.sort()

cutoff = 2
cycles = 5
#cgo_object = 1

cmd.load(target,target_name)

for i in range(len(file_list)):
	obj_name1 = extension.sub('',file_list[i])
	cmd.load(file_list[i],obj_name1)
	rms = cmd.super(obj_name1,target_name,cutoff=cutoff,cycles=cycles)
	print(str(obj_name1) + ' RMSD: ' + str(rms))
	print("saving pdb...")
	cmd.save(obj_name1 + "_aligned.pdb", obj_name1)
	cmd.delete(obj_name1)
