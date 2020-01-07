# ===============================================================================
#     This file is part of Jwalk.
#
#     Jwalk - A tool to calculate the solvent accessible surface distance (SASD)
#     between crosslinked residues.
#
#     Copyright 2016 Jwalk Inventor and Birkbeck College University of London.
#                          The Jwalk Inventor is: Josh Bullock
#
#
#     Jwalk is available under Public Licence.
#     This software is made available under GPL V3
#
#     Please cite your use of Jwalk in published work:
#
#     J.Bullock, J. Schwab, K. Thalassinos, M. Topf (2016)
#     The importance of non-accessible crosslinks and solvent accessible surface distance
#     in modelling proteins with restraints from crosslinking mass spectrometry.
#     Molecular and Cellular Proteomics (15) pp.2491-2500
#
# ===============================================================================

from Jwalk import PDBTools, GridTools, SurfaceTools, SASDTools
import math
import os
import sys
import argparse
from multiprocessing import cpu_count
import itertools

class ResidueClashException(Exception):
    """Raised when confusion about residue specification."""
    def __init__(self,msg):
        super().__init__(msg)

# default parameters

max_dist = 60
vox = 1
surface = False
surface_norm = False
surface_depth = False
bs3 = False
edc = False
xl_list = []
aa1 = "LYS"
aa2 = "LYS"
ncpus = cpu_count()
pdb_list = [i for i in os.listdir('./') if i.endswith('.pdb')]
bs3_combs = [i for i in itertools.combinations_with_replacement(['LYS','THR','SER','TYR'],2)]
edc_combs = [('LYS','ASP'),('LYS','GLU')]

amino_acids = {"LYS":"lysines",
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
                "TYR":"tyrosines"}

parser = argparse.ArgumentParser(description='JWALK: Calculate SASDs on your target PDB files')

parser.add_argument('-lys', action="store_true",
					help='calculate lysine crosslinks (default)')
parser.add_argument('-bs3', action="store_true",
					help='calculate lysine, tyrosine, threonine and serine crosslinks (for BS3 crosslinking)')
parser.add_argument('-edc', action="store_true",
					help='calculate lysine-aspartic/glutamic acid crosslinks (for EDC crosslinking)')
parser.add_argument('-xl_list', nargs=1,
					help='calculate crosslinks from input list')
parser.add_argument('-i', nargs=1,
					help='specify input pdb: -i <inputfile.pdb>')
parser.add_argument('-aa1', nargs= "*",
					help='specify start amino acid (three letter code e.g. LYS) - Multiple 3 letter AAs can be spearated by a comma. E.g. LYS,SER,THR,TYR')
parser.add_argument('-aa2', nargs= "*",
					help='specify end amino acid (three letter code e.g. LYS) - Multiple 3 letter AAs can be spearated by a comma. E.g. LYS,SER,THR,TYR')
parser.add_argument('-max_dist', nargs=1,
					help='specify maximum crosslink distance in Angstroms')
parser.add_argument('-vox', nargs=1,
					help='specify voxel size of grid')
parser.add_argument('-surface', action="store_true",
					help='use higher accuracy method to calculate solvent accessibility - requires Freesasa installation')
parser.add_argument('-ncpus', nargs=1,
					help='specify number of cpus to use')

args = parser.parse_args()

if args.lys:
	aa1 = "LYS"
	aa2 = "LYS"

elif (args.bs3 or args.edc) and (args.aa1 or args.aa2):
    raise ResidueClashException("Cannot use custom AA lists with predifined lists (bs3/edc).")

elif args.xl_list and (args.bs3 or args.edc or args.aa1 or args.aa2):
    raise ResidueClashException("Cannot use XL List with custom AA lists or predifined lists (bs3/edc).")

elif args.xl_list:
	xl_list = args.xl_list[0]

elif args.bs3 and args.edc:
    bs3 = True
    edc = True

elif args.bs3:
    bs3 = True

elif args.edc:
    edc = True

elif args.aa1:
    if args.aa2:
        aa1 = args.aa1[0].upper().split(',')
        if len(aa1) == 1:
            aa1 = aa1[0]
        aa2 = args.aa2[0].upper().split(',')
        if len(aa2) == 1:
            aa2 = aa2[0]
        # catch any dodgy typing
        if isinstance(aa1,str) and isinstance(aa2,str):
            if aa1 not in amino_acids or aa2 not in amino_acids:
                print("ERROR: Please type amino acid in three letter code format")
                sys.exit(2)
        else:
            err_test = []
            if isinstance(aa1,str):
                err_test.append(aa1)
            else:
                for a in aa1:
                    err_test.append(a)
            if isinstance(aa2,str):
                err_test.append(aa1)
            else:
                for a in aa2:
                    err_test.append(a)

            for a in err_test:
                if a not in amino_acids:
                    print("ERROR: Please type amino acid in three letter code format")
                    sys.exit(2)

    else:
        print("Please specify both aa1 AND aa2 if you want to use this option")
        sys.exit(2)

if args.max_dist:
	max_dist = int(args.max_dist[0])

if args.i:
	pdb_list = [args.i[0]]

if args.vox:
	vox = int(args.vox[0])

if args.surface:
	ACCESS_BIN = "freesasa"
	surface = True

if args.ncpus:
	ncpus = int(args.ncpus[0])

#######################################

def runJwalk(max_dist, vox, surface, surface_norm, surface_depth, xl_list, aa1, aa2, ncpus, pdb_list):
    """
        Execute Jwalk with processed command line options

			max_dist: maximum distance Jwalk will search
			vox: angstoms per voxel in grid
			surface: if True use higher resolution surface method
			xl_list: list of specific crosslinks to calculate
			aa1: starting residue type
			aa2: ending residues type
			ncpus: number of cpus to use
			pdb_list: default is all pdbs in directory, unless otherwise stated

    """
    for pdb in pdb_list:

        if bs3 == False and edc == False and isinstance(aa1, str) and isinstance(aa2, str):
            print("calculating crosslinks on", pdb)
            # load pdb into Jwalk
            structure_instance = PDBTools.read_PDB_file(pdb)
            # generate grid of voxel size (vox) that encapsulates pdb
            grid = GridTools.makeGrid(structure_instance, vox)

            # mark C-alpha positions on grid
            if xl_list: # if specific crosslinks need to be calculated
                crosslink_pairs, aa1_CA, aa2_CA = GridTools.mark_CAlphas_pairs(grid, structure_instance,  xl_list)
            else:
                crosslink_pairs = [] # na if searching every combination between residue types
                aa1_CA, aa2_CA = GridTools.markCAlphas(grid, structure_instance, aa1, aa2)

            if surface == True:
                # check more rigorously if residues are solvent accessible or not
                aa1_CA = SurfaceTools.check_solvent_accessibility_freesasa(pdb, aa1_CA, ACCESS_BIN, xl_list)
                if aa1 != aa2 or xl_list:
                    aa2_CA = SurfaceTools.check_solvent_accessibility_freesasa(pdb, aa2_CA, ACCESS_BIN, xl_list)
                else:
                    aa2_CA = aa1_CA.copy()

            dens_map = GridTools.generate_solvent_accessible_surface(grid, structure_instance, aa1_CA, aa2_CA)
            # identify which residues are on the surface
            aa1_voxels, remove_aa1 = GridTools.find_surface_voxels(aa1_CA, dens_map, surface, xl_list)
            aa2_voxels, remove_aa2 = GridTools.find_surface_voxels(aa2_CA, dens_map, surface, xl_list)

            crosslink_pairs = SurfaceTools.update_crosslink_pairs(crosslink_pairs, aa1_CA, aa2_CA, remove_aa1, remove_aa2)

            # calculate sasds
            sasds = SASDTools.parallel_BFS(aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA, crosslink_pairs,
                                           max_dist, vox, ncpus, xl_list)

            # remove duplicates
            sasds = GridTools.remove_duplicates(sasds)
            sasds = SASDTools.get_euclidean_distances(sasds, pdb, aa1, aa2)

            # output sasds to .pdb file and .txt file
            PDBTools.write_sasd_to_txt(sasds, pdb)
            PDBTools.write_sasd_to_pdb(dens_map, sasds, pdb)
            print(len(sasds), "SASDs calculated")
        else:
            print("calculating crosslinks on", pdb)
            # load pdb into Jwalk
            structure_instance = PDBTools.read_PDB_file(pdb)
            # generate grid of voxel size (vox) that encapsulates pdb
            grid = GridTools.makeGrid(structure_instance, vox)

            glob_sasds = {}

            if bs3 == True and edc == False:
                res_list = bs3_combs
            elif edc == True and bs3 == False:
                res_list = edc_combs
            elif bs3 == True and edc == True:
                res_list = bs3_combs + edc_combs
            elif bs3 == False and edc == False and isinstance(aa1,list) and isinstance(aa2,list):
                res_list =  itertools.product(aa1,aa2)
            else:
                res_list = False

            if res_list:
                for aa1,aa2 in res_list:
                    print(aa1,aa2)
                    # mark C-alpha positions on grid
                    if xl_list: # if specific crosslinks need to be calculated
                        crosslink_pairs, aa1_CA, aa2_CA = GridTools.mark_CAlphas_pairs(grid, structure_instance,  xl_list)
                    else:
                        crosslink_pairs = [] # na if searching every combination between residue types
                        aa1_CA, aa2_CA = GridTools.markCAlphas(grid, structure_instance, aa1, aa2)

                    if surface == True:
                        # check more rigorously if residues are solvent accessible or not
                        aa1_CA = SurfaceTools.check_solvent_accessibility_freesasa(pdb, aa1_CA, ACCESS_BIN, xl_list)
                        if aa1 != aa2 or xl_list:
                            aa2_CA = SurfaceTools.check_solvent_accessibility_freesasa(pdb, aa2_CA, ACCESS_BIN, xl_list)
                        else:
                            aa2_CA = aa1_CA.copy()

                    dens_map = GridTools.generate_solvent_accessible_surface(grid, structure_instance, aa1_CA, aa2_CA)
                    # identify which residues are on the surface
                    aa1_voxels, remove_aa1 = GridTools.find_surface_voxels(aa1_CA, dens_map, surface, xl_list)
                    aa2_voxels, remove_aa2 = GridTools.find_surface_voxels(aa2_CA, dens_map, surface, xl_list)

                    crosslink_pairs = SurfaceTools.update_crosslink_pairs(crosslink_pairs, aa1_CA, aa2_CA, remove_aa1, remove_aa2)

                    # calculate sasds
                    sasds = SASDTools.parallel_BFS(aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA, crosslink_pairs,
                                                   max_dist, vox, ncpus, xl_list)

                    # remove duplicates
                    sasds = GridTools.remove_duplicates(sasds)
                    sasds = SASDTools.get_euclidean_distances(sasds, pdb, aa1, aa2)

                    for k,v in sasds.items():
                        glob_sasds[k] = v

                # output sasds to .pdb file and .txt file
                PDBTools.write_sasd_to_txt(glob_sasds, pdb)
                PDBTools.write_sasd_to_pdb(dens_map, glob_sasds, pdb)
                print(len(glob_sasds), "SASDs calculated")

if __name__ == "__main__":
    runJwalk(max_dist, vox, surface, surface_norm, surface_depth, xl_list, aa1, aa2, ncpus, pdb_list)
