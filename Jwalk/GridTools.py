#===============================================================================
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
#===============================================================================

from numpy import array,  zeros, real,sqrt,exp
from collections import OrderedDict
import sys

class Map:
    """

    A class representing all information from a density map file.
    NOTE: Currently it can only read the CCP4/MRC  format.

    """

    def __init__(self, fullMap, origin, apix, filename, header=[]):
        """

        Read a map and its parameters in to Map class instance.

        *filename*
            name of map file.
        *origin*
            origin co-ordinates of the map (x_origin, y_origin, z_origin).
        *apix*
            grid spacing of map.
        *filename*
            filename of the Map instance

            NOTE: The *filename* 'build++copy' is reserved for copying of other Map class instances."""
        self.header = header
        self.origin = origin
        self.apix = apix
        self.filename = filename
        self.fullMap = fullMap

    def copy(self):
        """

        Return:
            copy of the Map.

        """
        copy = Map(self.fullMap.copy(), self.origin[:], self.apix, self.filename, self.header[:])
        return copy

    def box_size(self):
        """

        Return:
            size of the map array, in ZYX format.

        """
        return self.fullMap.shape

    def x_size(self):
        """

        Return:
            x size of the map array in x direction.

        """
        return self.fullMap.shape[2]

    def y_size(self):
        """

        Return:
            y size of the map array in y direction.

        """
        return self.fullMap.shape[1]

    def z_size(self):
        """

        Return:
            z size of the map array in z direction.

        """
        return self.fullMap.shape[0]

def makeGrid(struct, apix, resolution = 3, filename = "None"):
    """

    Returns protein grid.

    Arguments:

       *struct*
           Tempy structure instance
       *apix*
           angstroms per voxel

    """
    # Build empty template map based on the size of the protein and the resolution.
    extr = struct.get_extreme_values()
    edge = int(2*resolution/apix)+2
    x_size = int((extr[1]-extr[0])/apix)+edge
    y_size = int((extr[3]-extr[2])/apix)+edge
    z_size = int((extr[5]-extr[4])/apix)+edge

    # Origin calculated such that the centre of the map is the centre of mass of the protein.
    x_origin = (extr[1]+extr[0])/2-(apix*x_size/2.0)
    y_origin = (extr[3]+extr[2])/2-(apix*y_size/2.0)
    z_origin = (extr[5]+extr[4])/2-(apix*z_size/2.0)

    newMap = zeros((z_size, y_size, x_size))
    fullMap = Map(newMap, [x_origin, y_origin, z_origin], apix, filename)
    return fullMap

def mapGridPosition(densMap, atom):

    """

    Returns the index of the nearest pixel to an atom, and atom mass (4 values in list form).

    Arguments:

       *densMap*
           Map instance the atom is to be placed on.
       *atom*
           Atom instance.

    """
    origin = densMap.origin
    apix = densMap.apix
    box_size = densMap.box_size()
    x_pos = int(round((atom.x-origin[0])/apix,0))
    y_pos = int(round((atom.y-origin[1])/apix,0))
    z_pos = int(round((atom.z-origin[2])/apix,0))

    if((densMap.x_size() > x_pos >= 0) and (densMap.y_size() > y_pos >= 0) and (densMap.z_size() > z_pos >= 0)):
        return (x_pos, y_pos, z_pos, atom.mass)
    else:
        return 0

def markCAlphas(densMap, prot, aa1, aa2):
    """

    Returns ordered dictionaries containing {residue_number, chain, residue name : x, y, z}
    for both aa1 and aa2.

    Arguments:

       *densMap*
           Protein grid
       *prot*
           Tempy structure instance
        *aa1*
            Residue type 1
        *aa2*
            Residue type 2

    """

    aa1_CA = OrderedDict()
    aa2_CA = OrderedDict()

    for atom in prot.atomList:

            if atom.res == aa1:

                if atom.atom_name == 'CA':
                    pos = mapGridPosition(densMap, atom)
                    aa1_CA[atom.res_no,atom.chain,atom.res]=[pos[0],pos[1],pos[2]]

            if atom.res == aa2:
                if atom.atom_name == 'CA':
                    pos = mapGridPosition(densMap, atom)
                    aa2_CA[atom.res_no,atom.chain,atom.res]=[pos[0],pos[1],pos[2]]


    return aa1_CA, aa2_CA

def process_input_crosslinks(uv_xl):
    """Processes crosslink input .txt file and returns list of residues and crosslinked pairs"""

    aa1 = []
    aa2 = []
    crosslink_pairs = []
    #c = 0
    count = 0
    with open(uv_xl) as xl_in:
        for line in xl_in:
            count +=1
            col = line.split("|")
            try:
                chain1 = col[1].rstrip()
                chain1 = chain1.lstrip()
                chain2 = col[3].rstrip()
                chain2 = chain2.lstrip()
            except:
                print("ERROR: formatting error on line "+str(count)+" : "+line)
                exit(1)
            # if no chain is given
            if len(chain1) == 0:
                chain1 = " "
            if len(chain2) == 0:
                chain2 = " "

            aa1.append([int(col[0]),chain1])
            aa2.append([int(col[2]),chain2])

            crosslink_pairs.append([(int(col[0]),chain1),(int(col[2]),chain2)])
            '''
            if (int(col[0]),chain1,c) in xl_pair:
                c +=1
            xl_pair[int(col[0]),chain1,c] = [int(col[2]),chain2,c]
            '''

    return aa1, aa2 , crosslink_pairs

def mark_CAlphas_pairs(densMap, prot, uv_xl):

    """
    Processes input txt file. Checks each atom is in the structure and returns
    the pairs of crosslinks as well as the Calpha positions on the grid.

    Arguments:
        *densMap*
            grid that encompasses protein
        *prot*
            .pdb file
        *uv_xl*
             .txt input file

    """

    # process txt file
    aa1, aa2, crosslink_pairs = process_input_crosslinks(uv_xl)
    # crosslink_pairs is a list of [crosslinked aa1, aa2]

    aa1_CA = OrderedDict()
    aa2_CA = OrderedDict()

    atom_check = []

    for atom in prot.atomList:
        if [atom.res_no,atom.chain] in aa1:
            if atom.atom_name == "CA":
                pos = mapGridPosition(densMap, atom)
                aa1_CA[atom.res_no,atom.chain,atom.res]=[pos[0],pos[1],pos[2]]
                atom_check.append([atom.res_no,atom.chain])
        if [atom.res_no,atom.chain] in aa2:
            if atom.atom_name == "CA":
                pos = mapGridPosition(densMap, atom)
                aa2_CA[atom.res_no,atom.chain,atom.res]=[pos[0],pos[1],pos[2]]
                atom_check.append([atom.res_no,atom.chain])

    # check that all the residues listed are in the structure
    rem_x = []

    for x in aa1:
        if x not in atom_check:
            print("ERROR ! Residue",x[0],"-",x[1], "not in pdb structure - please check input files")
            rem_x.append((x[0], x[1]))
    for x in aa2:
        if x not in atom_check:
            print("ERROR ! Residue",x[0],"-",x[1], "not in pdb structure - please check input files")
            rem_x.append((x[0], x[1]))

    # remove crosslinks from crosslink_pairs if one or both residues are not in structure
    index_to_delete = []

    for i in range(len(crosslink_pairs)):
        [x1, x2] = crosslink_pairs[i]
        if x1 in rem_x:
            index_to_delete.append(i)
        elif x2 in rem_x:
            index_to_delete.append(i)

    crosslink_pairs_hold = []
    for i in range(len(crosslink_pairs)):
        if i not in index_to_delete:
            crosslink_pairs_hold.append(crosslink_pairs[i])

    # append residue name to crosslink_pairs_final

    aa_d = {}

    for a in aa1_CA:
        aa_d[a[0],a[1]]= a
    for a in aa2_CA:
        aa_d[a[0],a[1]]= a

    crosslink_pairs_final = []
    for x1, x2 in crosslink_pairs_hold:
        crosslink_pairs_final.append([aa_d[x1],aa_d[x2]])

    return crosslink_pairs_final, aa1_CA, aa2_CA

def generate_solvent_accessible_surface(densMap,prot,aa1_CA, aa2_CA):

    """

    Returns masked array which functions as solvent accessible surface

    Arguments:

       *densMap*
           Protein Grid
       *prot*
           Tempy structure instance
        *aa1_CA*
           voxel positions of each C_alpha atom of interest
        *aa2_CA*
           voxel positions of each C_alpha atom of interest

    """
    # store different radii for each atom and calculate the voxel spheres for each
    sphere = {}
    radius = {}

    C = 0.8

    radius['CA'] = 1.73 + C
    radius['S'] = 1.67 + C
    radius['N'] = 1.43 + C
    radius['OH'] = 1.30 + C

    for r in radius:

        sphere[r] = []
        rad = int(round(radius[r]/densMap.apix))

        for x in range(-rad,rad+1):
            for y in range(-rad,rad+1):
                for z in range(-rad,rad+1):
                    if (x**2 + y**2 + z**2) <= (rad**2):
                        sphere[r].append([x,y,z])

    backbone = ['N','CA','C','O']

    # generate solvent accessible surface

    for atom in prot.atomList:

        pos = mapGridPosition(densMap, atom)

        if pos:
            # don't place side chain atoms of residues of interest in the surface
            if ((atom.res_no,atom.chain,atom.res) in aa1_CA and atom.atom_name not in backbone) or (
            (atom.res_no,atom.chain,atom.res) in aa2_CA and atom.atom_name not in backbone):
                pass
            # for each atom, expand the corresponding voxel sphere around it to create solvent accessible surface
            else:
                if atom.atom_name[:1] == 'C':
                    for (x,y,z) in sphere['CA']:
                        if((densMap.x_size() > (pos[0]+x) >= 0) and (densMap.y_size() > (pos[1]+y) >= 0) and (densMap.z_size() > (pos[2]+z) >= 0)):
                            densMap.fullMap[pos[2]+z][pos[1]+y][pos[0]+x] += 1
                elif atom.atom_name[:1] == 'O':
                    for (x,y,z) in sphere['OH']:
                        if((densMap.x_size() > (pos[0]+x) >= 0) and (densMap.y_size() > (pos[1]+y) >= 0) and (densMap.z_size() > (pos[2]+z) >= 0)):
                            densMap.fullMap[pos[2]+z][pos[1]+y][pos[0]+x] += 1
                elif atom.atom_name[:1] == 'N':
                    for (x,y,z) in sphere['N']:
                        if((densMap.x_size() > (pos[0]+x) >= 0) and (densMap.y_size() > (pos[1]+y) >= 0) and (densMap.z_size() > (pos[2]+z) >= 0)):
                            densMap.fullMap[pos[2]+z][pos[1]+y][pos[0]+x] += 1
                elif atom.atom_name[:1] == 'S':
                    for (x,y,z) in sphere['S']:
                        if((densMap.x_size() > (pos[0]+x) >= 0) and (densMap.y_size() > (pos[1]+y) >= 0) and (densMap.z_size() > (pos[2]+z) >= 0)):
                            densMap.fullMap[pos[2]+z][pos[1]+y][pos[0]+x] += 1

    return densMap

def find_empty_space(res,sphere,densMap,CA):

    """

    Returns list of empty voxels in voxel sphere shell

    Arguments:

       *res*
           residue where search is happening around
       *sphere*
           voxel sphere shell to be expanded around CA voxel
        *densMap*
           Solvent accessible surface (masked array)
        *CA*
           Calpha voxels

    """

    starters = []
    (x,y,z) = CA[res]
    for (x_s,y_s,z_s) in sphere:
        if((densMap.x_size() > (x+x_s) >= 0) and (densMap.y_size() > (y+y_s) >= 0) and (densMap.z_size() > (z+z_s) >= 0)):
            if densMap.fullMap[z_s+z][y_s+y][x_s+x] <= 0:
                starters.append([x_s+x,y_s+y,z_s+z])
    return starters

def find_surface_voxels(aa1_CA, densMap, surface, xl_list = []):

    """

    Returns ordered dictionaries containing all possible staring voxels for each Calpha of
    interest. If Calpha is not solvent accessible then no starting voxels are returned.

    If xl_list flag True, then list of entries to be removed is also returned. empty list
    otherwise.

    Arguments:

       *aa1_CA*
           Calpha voxels for amino acid type 1
       *densMap*
           Solvent accessible surface (masked array)

    """

    # generate voxel spheres shells to progressively extend search for starting voxels
    # this is to keep starting voxels as close to CA as possible.

    sphere1 = []
    sphere2 = []
    sphere3 = []
    sphere4 = []
    sphere5 = []
    sphere6 = []

    C = 1.68

    radius = 1.73 + C

    radius4 = int(round(radius/densMap.apix))+1 # radius rounded up = 4 with apix = 1
    radius3 = radius4 -1
    radius2 = radius3 -1
    radius1 = radius2 -1
    radius5 = radius4 + 1
    radius6 = radius5 + 1

    for x in range(-radius1,radius1+1):
        for y in range(-radius1,radius1+1):
            for z in range(-radius1,radius1+1):
                if (x**2 + y**2 + z**2) <= (radius1**2):
                    sphere1.append([x,y,z])

    for x in range(-radius2,radius2+1):
        for y in range(-radius2,radius2+1):
            for z in range(-radius2,radius2+1):
                if (x**2 + y**2 + z**2) <= (radius2**2):
                    if ([x,y,z]) not in sphere1:
                        sphere2.append([x,y,z])

    for x in range(-radius3,radius3+1):
        for y in range(-radius3,radius3+1):
            for z in range(-radius3,radius3+1):
                if (x**2 + y**2 + z**2) <= (radius3**2):
                    if ([x,y,z]) not in sphere1 and ([x,y,z]) not in sphere2:
                        sphere3.append([x,y,z])

    for x in range(-radius4,radius4+1):
        for y in range(-radius4,radius4+1):
            for z in range(-radius4,radius4+1):
                if (x**2 + y**2 + z**2) <= (radius4**2):
                    if ([x,y,z]) not in sphere1 and ([x,y,z]) not in sphere2 and (
                    [x,y,z]) not in sphere3:
                        sphere4.append([x,y,z])

    for x in range(-radius5,radius5+1):
        for y in range(-radius5,radius5+1):
            for z in range(-radius5,radius5+1):
                if (x**2 + y**2 + z**2) <= (radius5**2):
                    if ([x,y,z]) not in sphere1 and ([x,y,z]) not in sphere2 and (
                    [x,y,z]) not in sphere3 and ([x,y,z]) not in sphere4:
                        sphere5.append([x,y,z])

    for x in range(-radius6,radius6+1):
        for y in range(-radius6,radius6+1):
            for z in range(-radius6,radius6+1):
                if (x**2 + y**2 + z**2) <= (radius6**2):
                    if ([x,y,z]) not in sphere1 and ([x,y,z]) not in sphere2 and (
                    [x,y,z]) not in sphere3 and ([x,y,z]) not in sphere4 and ([x,y,z]) not in sphere5:
                        sphere6.append([x,y,z])

    # iterate through sphere shells and append starting voxels to dictionary

    aa_1_start_voxels = OrderedDict()
    buried = []
    k_count = 0
    k_buried = 0

    for k in aa1_CA:
        k_count +=1
        aa_1_start_voxels[k] = find_empty_space(k,sphere1,densMap,aa1_CA)
        if aa_1_start_voxels[k] == []:
            aa_1_start_voxels[k] = find_empty_space(k,sphere2,densMap,aa1_CA)
        if aa_1_start_voxels[k] == []:
            aa_1_start_voxels[k] = find_empty_space(k,sphere3,densMap,aa1_CA)
        if aa_1_start_voxels[k] == []:
            aa_1_start_voxels[k] = find_empty_space(k,sphere4,densMap,aa1_CA)

        if surface == False:
            # if no voxels found after fourth shell then residue is not solvent exposed
            if aa_1_start_voxels[k] == []:
                buried.append(k)
                del aa_1_start_voxels[k]
                k_buried +=1
        else: # residues that have been determined solvent accessible by other means.
            if aa_1_start_voxels[k] == []:
                aa_1_start_voxels[k] = find_empty_space(k,sphere5,densMap,aa1_CA)
            if aa_1_start_voxels[k] == []:
                aa_1_start_voxels[k] = find_empty_space(k,sphere6,densMap,aa1_CA)

    rem_x = []

    if len(buried) > 0 and len(xl_list) > 0:
        print("ERROR - ",k_buried," buried residue(s) in xl_list:")
        for t in buried:
            print(str(t[0])+"-"+str(t[1])+"-"+str(t[2]))
            rem_x.append([t[0],t[1]])
        #sys.exit(2)

    return aa_1_start_voxels, rem_x

def remove_duplicates(sasds,xl_list=False):

    """

    Returns sasds with duplicates removed.

    Arguments:

       *sasds*
           dictionary of sasds

    """

    if xl_list == True:
        return sasds

    keep = {}
    keep_keys = []
    keep_sasds = {}

    for (start, end, distance) in sasds:

        # check if there is a duplicate
        if (start, end) not in keep:
            keep[start, end] = distance

        # if there is a duplicate check which has the shortest distance and keep
        elif (start, end) in keep:
            if (distance < keep[start, end]):
                keep[start, end] = distance

    # reform the dictionary key
    for (j,k) in keep:
        keep_keys.append((j,k,keep[j,k]))

    # filter out original dictionary
    for k in keep_keys:
        keep_sasds[k] = sasds[k]

    return keep_sasds
