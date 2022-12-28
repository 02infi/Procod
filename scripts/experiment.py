#!usr/bin/python
import numpy as np
from numpy import linalg as la 
import os
import math
import re 


def openfile(path): 
	info = []
	with open(path,"rU") as infile : 
		line_list = infile.readlines()
		for line in line_list:
			info.append(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",line))
	infile.close()
	return info


def complex():
	path_for_complex = "/home/gaurav/Gaurav/Project_1/Area/gepol/geopl-old/cpa/exp_4/complex.text"
	complexA = openfile(path_for_complex)
	return complexA


def moleculeAB():
	path_for_moleculeAB = "/home/gaurav/Gaurav/Project_1/Area/gepol/geopl-old/cpa/exp_4/molecule.text"
	molAB = openfile(path_for_moleculeAB)
	return molAB 


def normal_vector():
	vector_mol = (openfile("/home/gaurav/Gaurav/Project_1/Area/gepol/geopl-old/cpa/exp_4/complex.normal"))
	return vector_mol


def complex_pdb():
	info = []
	with open("/home/gaurav/Gaurav/Project_1/Area/gepol/geopl-old/cpa/exp_4/new.pdb","rU") as infile : 
		line_list = infile.readlines()
		for line in line_list: 
			info.append((line).split())
	infile.close()
	return info


def interface():
	data = []
	com = np.array(complex())
	mol = np.array(moleculeAB())
	row = np.array(complex()).shape[0]
	i = 0
	condition = True 

  	for eachrow in range(1,row):
  		if not (float (mol[eachrow,1]) == float (com[eachrow,1])):
  			if (float (com[eachrow,1]) < 5):
  				data.append(mol[eachrow,0])
  				if (int (mol[eachrow,0]) > i and condition):
  					i = int (mol[eachrow,0])
  				else : condition = False 
 	return data,i	


def penetration ():
	data,end_of_first_mol = interface()


 	
def normal_calculation():
	normal = []
	rectify = []
	resultant_molA = []
	resultant_molB = []
	mol = (normal_vector())
	for eachrow in range(1,len(mol)):
		if (mol[eachrow][1] == '  1'):
			rectify.append(mol[eachrow])
		else : continue
	
	data,i = interface()
	print data
	length = len(rectify)
	end = len(data)
	start_with_zero = 0
	condition = True

	for row in range(0,end):
		for lne in range(start_with_zero,length):
			if (int (data[row]) == int (rectify[lne][0])):
				normal = rectify[lne][-7:]
				normal.append(rectify[lne][0])

				if (int(rectify[lne][0]) <= i and condition) :
					resultant_molA.append(normal)
				else : 
					resultant_molB.append(normal) 
					condition = False 
				start_with_zero = lne
			else : continue  


	return np.array(resultant_molA),np.array(resultant_molB)


def corresponding_molecule(normal_A,normal_B):

	mol = np.array(moleculeAB())
	i = 0
	condition = True 

  	for eachrow in range(1,len(mol)):
  		if (int (mol[eachrow,0]) > i):
  			i = int (mol[eachrow,0])
  		else : break 

  	print i , len(mol) 

	distance_table = []
	position_vec_AB = []
	atom_A = np.array(normal_A[:,7],dtype = int)
	atom_B = np.array(normal_B[:,7],dtype = int)
	pdb_AB = complex_pdb()
	example = np.matrix(normal_B[:,0])
	A = np.array(normal_A[:,0:3],np.float)
	B = np.matrix(normal_B[:,0:3],np.float)
	normal_vec_A = np.matrix(normal_A[:,4:7],np.float)
	atom_data,endofmol_A = interface()
	magnitude_A = (np.apply_along_axis(np.linalg.norm, 1, np.array(normal_A[:,4:7],np.float)))
	pdb = np.array(pdb_AB,dtype = object)
	temp = atom_A[0]
	uniq = []
	uniq = np.array(uniq,dtype = int)

	for row in range(0,len(normal_A)):
		magnitude_distance = (np.apply_along_axis(np.linalg.norm,1,(np.cross(np.matrix(np.subtract(B,A[row,0:3])),np.matrix(normal_A[row,4:7],np.float)))))
		distance = np.array(np.divide(magnitude_distance,magnitude_A[row]))
		values = np.array(np.where(distance < 2))
		if not values.size: 
			continue
		Each_A = np.matrix(normal_A[row,0:3],np.float)
		Each_B = np.matrix(normal_B[values,0:3],np.float)
		position_vec_AB = np.subtract(Each_A,Each_B)	
		normalA = np.matrix(normal_A[row,4:7],np.float)
		normalB = np.matrix(normal_B[values,4:7],np.float)
		dot_product_AB = np.append(np.sum(np.multiply(normalA,-position_vec_AB),axis=1),np.sum(np.multiply(normalB,position_vec_AB),axis=1),axis=1)
		s = np.array(np.where(dot_product_AB < 0))[:][0]
		rep_el = np.array(s[np.diff(s) == 0],np.int)
		if (temp == atom_A[row]):
			uniq = np.concatenate((uniq,atom_B[values[0][rep_el]]))
		else : 
			penetrate = np.unique(uniq)
			if penetrate.size: 
				print temp,pdb[temp-1]
			temp = atom_A[row]
			uniq = []
			uniq = np.array(uniq,dtype = int)
			uniq = np.concatenate((uniq,atom_B[values[0][rep_el]]))			
	return 0 

def main():
	normal_A,normal_B = normal_calculation()
	print normal_B
	corresponding_molecule(normal_A,normal_B)

main()
