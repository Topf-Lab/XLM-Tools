#===============================================================================
#	 This file is part of Jwalk.
#
#	 Jwalk - A tool to calculate the solvent accessible surface distance (SASD)
#	 between crosslinked residues.
#
#	 Copyright 2016 Jwalk Inventor and Birkbeck College University of London.
#						  The Jwalk Inventor is: Josh Bullock
#
#
#	 Jwalk is available under Public Licence.
#	 This software is made available under GPL V3
#
#	 Please cite your use of Jwalk in published work:
#
#	 J.Bullock, J. Schwab, K. Thalassinos, M. Topf (2016)
#	 The importance of non-accessible crosslinks and solvent accessible surface distance
#	 in modelling proteins with restraints from crosslinking mass spectrometry.
#	 Molecular and Cellular Proteomics (15) pp.2491-2500
#
#===============================================================================

import sys
from math import cos, sin
import subprocess
import os

def expand_points(atom,sphere):
	""" Exapnd the unit sphere around specific x,y,z point and returns surface points"""
	points = {}
	atom.point_list = []
	CH2 = 1.68
	radius = {"N":1.43,
			  "O":1.30,
			  "C":1.68,
			  "S":1.67,
			  }

	r = radius[atom.atom_name[0]] + CH2
	(x,y,z) = (atom.x, atom.y, atom.z)

	for s in sphere:
		x1 = x + s[0]*r
		y1 = y + s[1]*r
		z1 = z + s[2]*r
		points[x1,y1,z1] = 0

	return points

def create_unit_sphere():
	""" Generates a unit sphere with 30 points on the surface """

	unit_sphere = []

	unit_sphere.append( [0.0, 1.0, 0.0] )
	unit_sphere.append( [0.0, -1.0, 0.0] )

	nstep = 5
	PI = 3.1415926536
	theta = PI/nstep
	arc = theta

	for istep in range(nstep):
		istep = istep + 1	# to change range from 0--9 to 1--10
		y1 = cos(istep*theta)
		r2 = sin(istep*theta)
		ndot2= 2*PI*r2/arc   # the circumference at that radius / proportion of pi
		if ndot2 == 0.0:
			continue
		theta2 = 2*PI/ndot2
		for idot in range(int(ndot2)):
			idot = idot + 1  # to change range from 0-- to 1--
			x2 = r2*cos(idot*theta2)
			z2 = r2*sin(idot*theta2)
			unit_sphere.append( [x2, y1, z2] )

	return unit_sphere

def check_solvent_accessibility(prot,aa1_CA,xl_list = False):

	'''

	Checks solvent accessibility of residues in aa1_CA
	Returns aa1_CA of solvent accessible residues.

	Arguments
		*prot*
			Tempy structure instance
		*aa1_CA*
			residues of interest

	'''
	# dictionary for text output
	string_dict = {"LYS":"lysines",
				"CYS":"cysteines",
				"ASP":"acidic residues",
				"GLU":"acidic residues",
				"VAL":"valines",
				"ILE":"isoleucines",
				"LEU":"leucines",
				"ARG":"arginines",
				"PRO":"prolines",
				"GLY":"glycines",
				"ALA":"alanines",
				"TRP":"tryptophans",
				"PHE":"phenylalanines",
				"SER":"serines",
				"GLN":"glutamines",
				"HIS":"histidines",
				"MET":"methionines",
				"THR":"threonines",
				"ASN":"asparagines",
				"TYR":"tyrosines"
				}

	radius = {"N":1.43,
			  "O":1.30,
			  "C":1.68,
			  "S":1.67,
			  }
	sd_res = False
	# create sphere of 30 points to expand around atoms
	sphere = create_unit_sphere()

	SA_res = {}
	CH2 = 1.68

	# this is not very efficient ... freesasa implementation to come
	for atom in prot.atomList:
		atom.clash = 1
		if (atom.res_no,atom.chain,atom.res) in aa1_CA:
			sd_res = atom.res
			# generate 30 points in unit spehere around atom
			points = expand_points(atom,sphere)

			for p in points:
				(x,y,z) = (p[0],p[1],p[2])
				# for every other atom check if points intersects with it
				for atom2 in prot.atomList:
					if atom2.res != atom.res or (atom2.res == atom.res and atom2.atom_name != atom.atom_name):
						r = radius[atom2.atom_name[0]] + CH2
						# need to transpose x,y,z
						(tx,ty,tz) = (x-atom2.x,y-atom2.y,z-atom2.z)
						# if the point lies within the sphere of that atom then it clashes
						if tx**2 + ty**2 + tz**2 <= r**2:
							points[p] = 1
							break
				# if any point on the sphere doesn't intersect with another then the atom is solvent accessible
				if points[p] == 0:
					atom.clash = 0
					break

			# if atom doesn't clash then residue information is kept in SA_res
			if atom.clash == 0:
				SA_res[atom.res_no,atom.chain,atom.res] = aa1_CA[atom.res_no,atom.chain,atom.res]

	# inform user on buried resiudes
	if xl_list:
		pass
	elif sd_res == "LYS":
		print("%d %s and 1 N-terminus of which %d are on the surface" % (len(aa1_CA)-1,string_dict[sd_res], len(SA_res)))
	else:
		print("%d %s of which %d are on the surface" % (len(aa1_CA),string_dict[sd_res], len(SA_res)))

	return SA_res

def update_crosslink_pairs(crosslink_pairs, aa1_CA, aa2_CA, remove_aa1, remove_aa2):

	'''Removes buried residues from crosslink_pairs'''

	buried_residues = []
	index_to_delete = []

	for i in range(len(crosslink_pairs)): # for each residue pair, check both are solvent accessible

		x1, x2 = crosslink_pairs[i]

		if x1 not in aa1_CA:
			index_to_delete.append(i)
			if x1 not in buried_residues:
				buried_residues.append(x1)
			if x2 not in aa2_CA and x2 not in buried_residues:
				buried_residues.append(x2)
		elif x2 not in aa2_CA:
			index_to_delete.append(i)
			if x2 not in buried_residues:
				buried_residues.append(x2)

		if [x1[0],x1[1]] in remove_aa1:
			index_to_delete.append(i)
			if x1 not in buried_residues:
				buried_residues.append(x1)
			if x2 in remove_aa2 and not x2 in buried_residues:
				buried_residues.append(x2)

		elif [x2[0],x2[1]] in remove_aa2:
			index_to_delete.append(i)
			if x2 not in buried_residues:
				buried_residues.append(x2)

	no_sasd_possible = []
	crosslink_pairs_final = []
	for i in range(len(crosslink_pairs)):
		if i not in index_to_delete:
			crosslink_pairs_final.append(crosslink_pairs[i])
		else:
			no_sasd_possible.append(crosslink_pairs[i])

	if len(no_sasd_possible) > 0:
		print("the following crosslinks cannot be calculated:")
		for s in no_sasd_possible:
			print("%s-%s-%s - %s-%s-%s" % (s[0][2],s[0][0],s[0][1],s[1][2],s[1][0],s[1][1]))

	return crosslink_pairs_final

def check_solvent_accessibility_freesasa(prot,aa1_CA,freesasa_source = "freesasa", xl_list = False):
	try:
		os.mkdir('rsa_files')
	except:
		pass
	# dictionary for text output
	string_dict = {"LYS":"lysines",
				"CYS":"cysteines",
				"ASP":"aspartates",
				"GLU":"glutamates",
				"VAL":"valines",
				"ILE":"isoleucines",
				"LEU":"leucines",
				"ARG":"arginines",
				"PRO":"prolines",
				"GLY":"glycines",
				"ALA":"alanines",
				"TRP":"tryptophans",
				"PHE":"phenylalanines",
				"SER":"serines",
				"GLN":"glutamines",
				"HIS":"histidines",
				"MET":"methionines",
				"THR":"threonines",
				"ASN":"asparagines",
				"TYR":"tyrosines"
				}

	cmd = [freesasa_source,"--rsa","-o","rsa_files/%s.rsa" % prot[:-4],"--radii=naccess",'--no-log','--no-warnings','--error-file=fs_log.txt',prot,'>','fs_log.txt','&']

	ext_freesasa = subprocess.Popen(cmd)
	ext_freesasa.communicate()

	try:
		test_open = open("rsa_files/"+prot[:-4]+".rsa")
		test_open.close
	except:
		print("No .rsa file created. Please check you have Freesasa installed. Go to 'http://freesasa.github.io/' ")
		exit(1)

	#try:
	SA_res = {}
	with open("rsa_files/"+prot[:-4]+".rsa") as rsa:
		for line in rsa:
			if line.startswith("RES") and not line.startswith('RESULTS'):
				col = line.split()
				RES = line[4:8].strip()
				CHAIN = line[8:9].strip()
				if len(CHAIN) == 0:
					CHAIN = " "
				NUM = line[9:13].strip()
				REL = line[23:29].strip()
				if float(REL) > 7.0:

					res = RES + CHAIN + NUM
					SA_res[(int(NUM),CHAIN,RES)] = True
	#except:
	#	print "error"
	SA_res_dict = {}

	for s in aa1_CA:
		if (s[0],s[1],s[2]) in SA_res:
			SA_res_dict[(s[0],s[1],s[2])] = aa1_CA[(s[0],s[1],s[2])]
			sd_res = s[2]
		else:
			print("Residue %d-%s-%s is buried" % (s[0],s[1],s[2]))
			sd_res = s[2]

	# inform user on buried resiudes
	if xl_list:
		pass
	elif sd_res == "LYS":
		print("%d %s and 1 N-terminus of which %d are on the surface" % (len(aa1_CA)-1,string_dict[sd_res], len(SA_res_dict)))
	else:
		print("%d %s of which %d are on the surface" % (len(aa1_CA),string_dict[sd_res], len(SA_res_dict)))

	return SA_res_dict

def check_solvent_accessibility_norm_freesasa(prot,aa1_CA,freesasa_source = "freesasa", xl_list = False):
	try:
		os.mkdir('rsa_files')
	except:
		pass

	x = -1
	# dictionary for text output
	string_dict = {"LYS":"lysines",
				"CYS":"cysteines",
				"ASP":"aspartates",
				"GLU":"glutamates",
				"VAL":"valines",
				"ILE":"isoleucines",
				"LEU":"leucines",
				"ARG":"arginines",
				"PRO":"prolines",
				"GLY":"glycines",
				"ALA":"alanines",
				"TRP":"tryptophans",
				"PHE":"phenylalanines",
				"SER":"serines",
				"GLN":"glutamines",
				"HIS":"histidines",
				"MET":"methionines",
				"THR":"threonines",
				"ASN":"asparagines",
				"TYR":"tyrosines"
				}


	cmd = [freesasa_source,"--rsa","-o","rsa_files/%s.rsa" % prot[:-4],"--radii=naccess",'--no-log','--no-warnings','--error-file=fs_log.txt',prot,'>','fs_log.txt','&']
	#cmd = [freesasa_source,"--format=rsa","--output=%s.rsa" % prot[:-4],"--radii=naccess",prot]

	ext_freesasa = subprocess.Popen(cmd)
	ext_freesasa.communicate()

	try:
		test_open = open("rsa_files/"+prot[:-4]+".rsa")
		test_open.close
	except:
		print("No .rsa file created. Please check you have Freesasa installed. Go to 'http://freesasa.github.io/' ")
		exit(1)

	#try:
	SA_res = {}
	with open("rsa_files/"+prot[:-4]+".rsa") as rsa:
		g = rsa.read().splitlines()
		rels = []
		for line in g:
			if line.startswith("RES") and not line.startswith('RESULTS'):
				rels.append(float(line[23:29].strip()))
		norm_rels = [(float(i)-min(rels))/(max(rels)-min(rels)) for i in rels]

		for line in g:
			if line.startswith("RES") and not line.startswith('RESULTS'):
				x+=1
				col = line.split()
				RES = line[4:8].strip()
				CHAIN = line[8:9].strip()
				if len(CHAIN) == 0:
					## Changed MS 150318 ##
					CHAIN = " "
				NUM = line[9:13].strip()
				REL = norm_rels[x]
				if float(REL) > float(0.01):

					res = RES + CHAIN + NUM
					SA_res[(int(NUM),CHAIN,RES)] = True
	#except:
	#	print "error"
	SA_res_dict = {}

	for s in aa1_CA:
		if (s[0],s[1],s[2]) in SA_res:
			SA_res_dict[(s[0],s[1],s[2])] = aa1_CA[(s[0],s[1],s[2])]
			sd_res = s[2]
		else:
			print("Residue %d-%s-%s is buried" % (s[0],s[1],s[2]))
			sd_res = s[2]

	# inform user on buried resiudes
	if xl_list:
		pass
	elif sd_res == "LYS":
		print("%d %s and 1 N-terminus of which %d are on the surface" % (len(aa1_CA)-1,string_dict[sd_res], len(SA_res_dict)))
	else:
		print("%d %s of which %d are on the surface" % (len(aa1_CA),string_dict[sd_res], len(SA_res_dict)))

	return SA_res_dict
