#!/usr/bin/python

import sys,re
import numpy as np
from numpy import linalg as la 
import os
import math
import shlex,subprocess
import codecs

def interface(complexAB,moleculeAB):
	data = []
	com = np.array(complexAB)
	mol = np.array(moleculeAB)
	row = np.array(complexAB).shape[0]
	i = 0
	j = 0
	condition_i = 1 
	condition_j = 1

  	for eachrow in range(1,row):
  		if not (float (mol[eachrow,1]) == float (com[eachrow,1])):
  			if (float (com[eachrow,1])):
  				data.append(mol[eachrow,0:6])
  				if ((int (mol[eachrow,0]) > i) & (condition_i)):
  					
  					i = int (mol[eachrow,0])
  				else : condition_i = 0 

  		if ((int(mol[eachrow,0]) > j) & (condition_j)):
  			j = int (mol[eachrow,0])
  		else : condition_j = 0

 	return np.array(data),i,j	
 	
def normal_calculation(data,i,j,normal_vector):
	normal = []
	center_A = []
	center_B = []
	rectify = []
	resultant_molA = []
	resultant_molB = []
	mol = normal_vector
	for eachrow in range(1,len(mol)):
		if (int (mol[eachrow][1]) == 1):
			rectify.append(mol[eachrow])
		else : continue
	
	length = len(rectify)
	end = len(data)
	start_with_zero = 0
	Switch = 1 
	Switch_A = 1

	for row in range(0,end):
		for lne in range(start_with_zero,length):

			if (int (data[row,0]) == int (rectify[lne][0])):
				Switch = 0
				normal = rectify[lne][-7:]
				normal.append(rectify[lne][0])

				if (Switch_A == 1):
					resultant_molA.append(normal)
					center_A.append(data[row])
				else :
					resultant_molB.append(normal)
					center_B.append(data[row])

			elif (Switch == 0):
				if (int (data[row,0]) == i):
					Switch_A = 0
				Switch = 1				
				start_with_zero = lne
				break

			else : continue 

	return np.array(resultant_molA),np.array(resultant_molB),np.array(center_A),np.array(center_B),j



def corresponding_molecule(normal_A,normal_B,center_A,center_B,k,pdb,Chain):
	atoms = []
	penetrating_res = []
	position_vec_AB = []
	atom_A = np.array(normal_A[:,7],dtype = int)
	atom_B = np.array(normal_B[:,7],dtype = int)
	pdb_AB = pdb
	example = np.matrix(normal_B[:,0])
	#A = np.array(center_A[:,2:5],np.float)
	#B = np.array(center_B[:,2:5],np.float)
	A = np.array(normal_A[:,0:3],np.float)
	B = np.array(normal_B[:,0:3],np.float)
	magnitude_A = (np.apply_along_axis(np.linalg.norm, 1, np.array(normal_A[:,4:7],np.float)))
	pdb = np.array(pdb_AB,dtype = 'string')
	temp = atom_A[0]
	uniq = []
	uniq = np.array(uniq,dtype = int)
	Area = 0
	condition = 1
	ligand_atoms = []
	ligand_atoms = np.array(ligand_atoms,dtype = int)
	Triangle = []

	for row in range(0,len(normal_A)):
		magnitude_distance = (np.apply_along_axis(np.linalg.norm,1,(np.cross(np.matrix(np.subtract(B,A[row,0:3])),np.matrix(normal_A[row,4:7],np.float)))))
		distance = np.array(np.divide(magnitude_distance,magnitude_A[row]))
		values = np.array(np.where(distance < 0.3))
		if not values.size: 
			continue
		Each_A = np.matrix(A[row,0:3],np.float)
		Each_B = np.matrix(B[values,0:3],np.float)
		position_vec_AB = np.subtract(Each_A,Each_B)
		normal_AB = np.apply_along_axis(np.linalg.norm,1,position_vec_AB)
		posit_unit_vec_AB = np.divide(position_vec_AB,normal_AB.reshape(len(normal_AB),1))
		normalA = np.matrix(normal_A[row,4:7],np.float)
		normalB = np.matrix(normal_B[values,4:7],np.float)
		dot_product_AB = np.append(np.sum(np.multiply(normalA,-position_vec_AB),axis=1),np.sum(np.multiply(normalB,position_vec_AB),axis=1),axis=1)
		s = np.array(np.where(dot_product_AB < 0))[:][0] 
		rep_el = np.array(s[np.diff(s) == 0],np.int)
		column = values[0][rep_el]

		if rep_el.size: 
			atoms.append(int(normal_A[row,7]))
			#print atom_A[row],atom_B[values[0][rep_el]]
			value = np.where(np.add(np.matrix(center_A[row,5],float),np.matrix(center_B[column,5],float)) > np.apply_along_axis(np.linalg.norm,1,np.subtract(np.array(center_A[row,2:5],np.float),np.array(center_B[column,2:5],np.float))))[:][1]

			if value.size:	
				Triangle.append(normal_A[row])
				if (condition == 1):
					temp = atom_A[row]
					ligand_atoms = np.concatenate((ligand_atoms,temp.ravel())) 
				#print pdb[atom_A[row]-1],pdb[atom_B[column[value]]+k-1]
				#print atom_A[row],atom_B[column[value]]
					
				if (temp == atom_A[row]):
					uniq = np.concatenate((uniq,atom_B[column[value]].ravel()))
					condition = 0
					#print temp,uniq
				else :
					#protein_atoms = np.concatenate((protein_atoms,np.unique(uniq)))
					condition = 1 
					uniq = []
					uniq = np.array(uniq,dtype = int)
					row = row-1
	#protein_atoms = np.concatenate((protein_atoms,np.unique(uniq)))	
	infoA = []
	infoB = []

	for atom in range(0,len(ligand_atoms)):
		for cenA in range(0,len(center_A)):
			if(int (ligand_atoms[atom]) == int(center_A[cenA,0])):
				infoA.append(center_A[cenA])
				Area = Area + float(center_A[cenA,1])
				break

	#	for cenB in range(0,len(center_B)): 	
	#		if(int (protein_atoms[atom]) == int(center_B[cenB,0])):
	#			infoB.append(center_B[cenB])
	#			break
	Proteins = pdb[ligand_atoms-1]
	#print pdb[protein_atoms+k-1]
	infoA = np.array(infoA)
	#distances = np.add(np.matrix(infoA[:,5],float),np.matrix(infoB[:,5],float)) - np.apply_along_axis(np.linalg.norm,1,np.subtract(np.array(infoA[:,2:5],np.float),np.array(infoB[:,2:5],np.float)))
	#depth = np.amax(distances)
	#penetrating_volume(infoA,infoB,distances)
	#print "Penetrating Depth = %f Angstroms" %depth
	return Triangle,Proteins,Area 

def intermediate(normal_A,normal_B,center_A,center_B,j,pdb,Chain):
	if (normal_A.size & normal_B.size):
		Triangle,Proteins,Area = corresponding_molecule(normal_A,normal_B,center_A,center_B,j,pdb,Chain)
		if not len(Triangle):
			print "The molecules is not penetrating and its ready for calculating the respective potentials" 
	else: print "The molecules is not penetrating and its ready for calculating the respective potentials"
	return Triangle,Proteins,Area



def gepol(l_atom,l_residue):
	atom_data = []
	normal_data = []
	with open("data.text","w") as myfile:
		start = 1.0 
		gepol_input  = "input.text"
		gepol_output = "output.text"  
		myfile.write("%d\n%d\n%d\n%s\n%s\n%d"%(l_atom,start,l_residue,gepol_input,gepol_output,start)) 
	myfile.close()
	os.system("./pdbgepol <data.text")
	os.system("./gepol <output.text > gepol.text")
	os.remove("gepol.text")
	os.remove("output.text")
	os.remove("data.text")
	os.remove("fort.15")
	os.remove("input.text")
	with open("fort.7","r") as myfile:
		for _ in xrange(8):
			next(myfile)
		for line in myfile:
			atom_data.append(line.strip())
		myfile.close()
	os.remove("fort.7")
	with open("fort.8","r") as myfile:
		for _ in xrange(6):
			next(myfile)
		for line in myfile:
			normal_data.append(line.strip())
		myfile.close()
	os.remove("fort.8")
	return (np.array(atom_data),np.array(normal_data)) 


def formating(info):
	data = [] 
	length = len(info)
	atoms = 0
	resi_id = info[0][22:26]
	res_number = 1
	with open("input.text","w") as myfile:
		myfile.write("REMARK PCI\n") 
		for id in range(0,length):
			if (info[id][0:4] == "ATOM"):
				if(resi_id == info[id][22:26]):
					number = int(res_number)  
				else :
					res_number = res_number + 1
					number = int(res_number)
					resi_id = info[id][22:26] 

				myfile.write("%-9s %-12s %-6d %7s %7s %7s\n"%(info[id][0:6],info[id][12:17],number,info[id][30:38],info[id][38:46],info[id][46:54]))
				line = ("%6s%5s %4s%c%3s %c%4s%c   %8s%8s%8s"%(info[id][0:6],id+1,info[id][12:16],info[id][16],info[id][17:20],info[id][21],number,info[id][26],info[id][30:38],info[id][38:46],info[id][46:54]))
				data.append(line)
				atoms = atoms + 1
			else : continue 
		myfile.close()
	atom_data,normal_data = gepol(atoms,number)
	return (np.array(atom_data),np.array(normal_data),np.array(data)) 


def openfile(path): 
	info = []
	with open(path, "r") as fp:
		data = fp.readlines()
		for line in data :
			info.append(line.rstrip())
	fp.close()
	return np.array(info)

def statments(A,B,Chain):
	info_C = []
	info_M = []
	normal = []
	pdb =[]
	atom_A,normal_A,pdb_A  = formating(openfile(A))
	atom_B,normal_B,pdb_B = formating(openfile(B))
	molAB = np.concatenate((atom_A,atom_B))
	normal_vector = np.concatenate((normal_A,normal_B))
	pdbAB = np.concatenate((pdb_A,pdb_B))
	atom_AB,normal_AB,pdb_AB = formating(pdbAB)
	for index in range(0,len(atom_AB)):
		info_C.append(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",atom_AB[index]))
	complexAB = np.array(info_C)[np.where(np.array(info_C)[:,8] == "      0")[:][0]][:]
	for index in range(0,len(molAB)):
		info_M.append(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",molAB[index]))
	moleculeAB = np.array(info_M)[np.where(np.array(info_M)[:,8] == "      0")[:][0]][:]
	for index in range(0,len(normal_vector)):
		normal.append(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",normal_vector[index]))
	normalAB = normal
	for line in pdb_AB:
		pdb.append(line.strip())
	combine_pdb = pdb
	data,i,j = interface(complexAB,moleculeAB)
	normal_A,normal_B,center_A,center_B,j = normal_calculation(data,i,j,normalAB)
	Triangle,Proteins,Area = intermediate(normal_A,normal_B,center_A,center_B,j,pdbAB,Chain)
	#intermediate(normal_B,normal_A,center_B,center_A,0)
	return Triangle,Proteins,Area

def main(argv):
    if len(sys.argv)!= 3 :
    	print "Usage : python collision_detection.py A.pdb B.pdb"
    	sys.exit(1)
    else :
    	A,B = argv
    	Triangle_A , Proteins_A , Area_A = statments(A,B,1)
    	Triangle_B , Proteins_B , Area_B = statments(B,A,2)
    	Triangle = np.concatenate((Triangle_A,Triangle_B))
      	with open('Triangle.text', 'w') as f:
      		f.write(np.array2string(Triangle, separator=', '))
      	f.close()

    	Proteins = np.concatenate((Proteins_A,Proteins_B))
    	with open("Collision.log","w") as f:
    		f.write(np.array2string(Proteins, separator=', '))
    		f.write("\n")
		f.close()

		A = np.array(Area_A[:,5])
		B = np.array(Area_B[:,5])
		with open("Collision.log","w") as f:
			f.write("Primary Chain :")
			f.write(A)
			f.write("The Acessibility Surface Area of Penetrating Primary A : %f"%Area_A)
    		f.write("\n")
    		f.write("The Acessibility Surface Area of Penetrating Secondary A : %f"%Area_B)
		f.close




if __name__ == "__main__":
   main(sys.argv[1:])
