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

from multiprocessing import Pool, Process, freeze_support
import itertools
import sys
import math
import os

def calculate_specific_SASD(single_crosslink, aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA,
                            max_dist, vox):

    '''

    Breadth First Search of grid. For general info on algorithm see:
    https://en.wikipedia.org/wiki/Breadth-first_search

    Returns dictionary containing solvent accessible surface distances between specific starting res
    and ending res.

    {start res, end res, length in angstroms : voxel path of sasd}

    Arguments:

       *single_crosslink*
           start and end residue.
           start is key of aa1_voxels. aa1_voxels[start_residue] = all the starting voxels for that
           residue
       *aa1_voxels*
           dictionary containing starting voxels {start_residue : starting voxels}
       *aa2_voxels*
           dictionary containing ending voxels {end_residue : ending voxels}
       *dens_map*
           grid with solvent accessible surface (masked array)
       *aa1_CA*
           dictionary containing voxel of C-alpha
       *aa2_CA*
           dictionary containing voxel of C-alpha
       *max_dist*
           maximum distance BFS will search until
       *vox*
           number of angstoms per voxel

    '''

    start_residue = single_crosslink[0]
    end_residue = single_crosslink[1]

    specific_xl = {}

    comb = [[1, 0, 0],
            [-1, 0, 0],
            [0, -1, 0],
            [0, 1, 0],
            [0, 0, -1],
            [0, 0, 1],
            [1, 0, 1],
            [-1, 0, 1],
            [0, 1, 1],
            [0, -1, 1],
            [1, -1, 0],
            [-1, -1, 0],
            [1, 1, 0],
            [-1, 1, 0],
            [1, 0, -1],
            [-1, 0, -1],
            [0, 1, -1],
            [0, -1, -1],
            [1, 1, 1],
            [1, -1, 1],
            [-1, 1, 1],
            [-1, -1, 1],
            [1, 1, -1],
            [1, -1, -1],
            [-1, 1, -1],
            [-1, -1, -1]]

    # distance of diagonal steps
    diag1 = (math.sqrt((vox ** 2) * 2)) # 2d diagonal
    diag2 = (math.sqrt((vox ** 2) * 3)) # 3d diagonal

    queue = [] # voxels in queue for searching
    end_voxels = [] # list of voxels to find path to
    visited = {}  # list works as all the coordinates that have been visited - dictionary gives the path to said coordinate from startpoint
    distance = {} # keeps distance from starting voxel for each other voxel

    # place starting voxels into queue and initialise visited and distance
    for j in aa1_voxels[start_residue]:
        queue.append([j[0], j[1], j[2]])
        visited[j[0], j[1], j[2]] = [[j[0], j[1], j[2]]]
        distance[j[0], j[1], j[2]] = 0

    while queue:
        x_n, y_n, z_n = queue.pop(0)
        if distance[x_n, y_n, z_n] <= max_dist:
            for c in comb:
                x_temp = x_n + c[0]
                y_temp = y_n + c[1]
                z_temp = z_n + c[2]
                if (x_temp, y_temp, z_temp) not in visited:
                    if ((0 <= x_temp < dens_map.x_size()) and (0 <= y_temp < dens_map.y_size()) and (
                            0 <= z_temp < dens_map.z_size())):
                        temp_list = visited[x_n, y_n, z_n][:]
                        temp_list.append([x_temp, y_temp, z_temp])
                        visited[x_temp, y_temp, z_temp] = temp_list  # updated visited list

                        if dens_map.fullMap[z_temp][y_temp][x_temp] <= 0:  # if the voxel is in empty space
                            queue.append(([x_temp, y_temp, z_temp]))
                        # calculate the distance
                        diff_x = x_temp - x_n
                        diff_y = y_temp - y_n
                        diff_z = z_temp - z_n
                        if diff_x != 0 and diff_y != 0 and diff_z != 0:
                            distance[x_temp, y_temp, z_temp] = distance[x_n, y_n, z_n] + diag2
                        elif diff_x != 0 and diff_y != 0:
                            distance[x_temp, y_temp, z_temp] = distance[x_n, y_n, z_n] + diag1
                        elif diff_x != 0 and diff_z != 0:
                            distance[x_temp, y_temp, z_temp] = distance[x_n, y_n, z_n] + diag1
                        elif diff_y != 0 and diff_z != 0:
                            distance[x_temp, y_temp, z_temp] = distance[x_n, y_n, z_n] + diag1
                        else:
                            distance[x_temp, y_temp, z_temp] = distance[x_n, y_n, z_n] + vox

    # now we have a full set of paths into empty space starting from start_residue
    # all stored in visited. Now need to extract paths to specific residue
    shortest_distance = 9999
    all_distances = {}

    for j in aa2_voxels[end_residue]:

        (x, y, z) = j

        if (x, y, z) in visited:

            visited[(x, y, z)].insert(0, aa1_CA[start_residue]) # add aa1 CA voxel to path
            visited[(x, y, z)].append(aa2_CA[end_residue]) # add aa2 CA voxel to end of path

            # add the distance between starting/ending residue CA voxel and start/end voxel in path
            for i in [1, len(visited[(x, y, z)]) - 1]:
                (x_1, y_1, z_1) = visited[(x, y, z)][i - 1]
                (x_2, y_2, z_2) = visited[(x, y, z)][i]
                distance[(x, y, z)] += math.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2)

            all_distances[distance[(x, y, z)]] = visited[(x, y, z)]  # linking distance:path

            # keep record of shortest distance
            if shortest_distance > distance[(x, y, z)]:
                shortest_distance = distance[(x, y, z)]

            # now adding shortest xl to the final list

            if shortest_distance != 9999:
                # this is just to order the dict so that chain goes alphabetically
                specific_xl[start_residue, end_residue, shortest_distance] = all_distances[
                    shortest_distance]  # start lys, end lys, length of xl = path of xl

    return specific_xl


def calculate_SASDs(start_residue, aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA,
                    max_dist, vox):

    """

    Breadth First Search of grid. For general info on algorithm see:
    https://en.wikipedia.org/wiki/Breadth-first_search

    Returns dictionary containing solvent accessible surface distances between starting res
    and all possible ending res.

    {start res, end res, length in angstroms : voxel path of sasd}

    Arguments:

       *start_residue*
           key of aa1_voxels. aa1_voxels[start_residue] = all the starting voxels for that
           residue
       *aa1_voxels*
           dictionary containing starting voxels {start_residue : starting voxels}
       *aa2_voxels*
           dictionary containing ending voxels {end_residue : ending voxels}
       *dens_map*
           grid with solvent accessible surface (masked array)
       *aa1_CA*
           dictionary containing voxel of C-alpha
       *aa2_CA*
           dictionary containing voxel of C-alpha
       *max_dist*
           maximum distance BFS will search until
       *vox*
           number of angstoms per voxel


    """

    sasds = {}

    # order of voxels to search - by having diagonals last ensures shortest path is returned
    comb = [[1, 0, 0],
            [-1, 0, 0],
            [0, -1, 0],
            [0, 1, 0],
            [0, 0, -1],
            [0, 0, 1],
            [1, 0, 1],
            [-1, 0, 1],
            [0, 1, 1],
            [0, -1, 1],
            [1, -1, 0],
            [-1, -1, 0],
            [1, 1, 0],
            [-1, 1, 0],
            [1, 0, -1],
            [-1, 0, -1],
            [0, 1, -1],
            [0, -1, -1],
            [1, 1, 1],
            [1, -1, 1],
            [-1, 1, 1],
            [-1, -1, 1],
            [1, 1, -1],
            [1, -1, -1],
            [-1, 1, -1],
            [-1, -1, -1]]

    # distance of diagonal steps
    diag1 = (math.sqrt((vox ** 2) * 2)) # 2d diagonal
    diag2 = (math.sqrt((vox ** 2) * 3)) # 3d diagonal

    queue = [] # voxels in queue for searching
    visited = {}  # list works as all the coordinates that have been visited - dictionary gives the path to said coordinate from startpoint
    distance = {} # keeps distance from starting voxel for each other voxel

    # place starting voxels into queue and initialise visited and distance
    for j in aa1_voxels[start_residue]:
        queue.append([j[0], j[1], j[2]])
        visited[j[0], j[1], j[2]] = [[j[0], j[1], j[2]]]
        distance[j[0], j[1], j[2]] = 0

    # grid is searched until queue is empty
    while queue:
        x_n, y_n, z_n = queue.pop(0) # take first voxel in queue
        if distance[x_n, y_n, z_n] <= max_dist:
            for c in comb: # expand in all directions from voxel - in order of comb.
                x_temp = x_n + c[0]
                y_temp = y_n + c[1]
                z_temp = z_n + c[2]
                # check voxel hasn't already been searched
                if (x_temp, y_temp, z_temp) not in visited:
                    # check that voxel is within bounds of the grid
                    if ((0 <= x_temp < dens_map.x_size()) and (0 <= y_temp < dens_map.y_size()) and (
                            0 <= z_temp < dens_map.z_size())):
                        # add path to this voxel to visited
                        temp_list = visited[x_n, y_n, z_n][:]
                        temp_list.append([x_temp, y_temp, z_temp])
                        visited[x_temp, y_temp, z_temp] = temp_list

                        if dens_map.fullMap[z_temp][y_temp][x_temp] <= 0:  # if the voxel is in empty space
                            queue.append(([x_temp, y_temp, z_temp])) # add to queue for later searching

                        # calculate the distance to voxel from start voxel
                        diff_x = x_temp - x_n
                        diff_y = y_temp - y_n
                        diff_z = z_temp - z_n
                        if diff_x != 0 and diff_y != 0 and diff_z != 0:
                            distance[x_temp, y_temp, z_temp] = distance[x_n, y_n, z_n] + diag2
                        elif diff_x != 0 and diff_y != 0:
                            distance[x_temp, y_temp, z_temp] = distance[x_n, y_n, z_n] + diag1
                        elif diff_x != 0 and diff_z != 0:
                            distance[x_temp, y_temp, z_temp] = distance[x_n, y_n, z_n] + diag1
                        elif diff_y != 0 and diff_z != 0:
                            distance[x_temp, y_temp, z_temp] = distance[x_n, y_n, z_n] + diag1
                        else:
                            distance[x_temp, y_temp, z_temp] = distance[x_n, y_n, z_n] + vox

    # now we have a full set of paths into empty space starting from start_residue
    # all stored in visited. Now need to extract paths to specific residues
    for end_residue in aa2_voxels:
        if start_residue != end_residue:
            shortest_distance = 9999
            all_distances = {}

            # cycling through possible end coords of end_residue to get shortest sasd
            for j in aa2_voxels[end_residue]:

                (x, y, z) = j

                if (x, y, z) in visited:

                    visited[(x, y, z)].insert(0, aa1_CA[start_residue]) # add aa1 CA voxel to path
                    visited[(x, y, z)].append(aa2_CA[end_residue]) # add aa2 CA voxel to end of path

                    # add the distance between starting/ending residue CA voxel and start/end voxel in path
                    for i in [1, len(visited[(x, y, z)]) - 1]:
                        (x_1, y_1, z_1) = visited[(x, y, z)][i - 1]
                        (x_2, y_2, z_2) = visited[(x, y, z)][i]
                        distance[(x, y, z)] += math.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2)

                    all_distances[distance[(x, y, z)]] = visited[(x, y, z)]  # linking distance:path

                    # keep record of shortest distance
                    if shortest_distance > distance[(x, y, z)]:
                        shortest_distance = distance[(x, y, z)]

            # add shortest distance sasd to output dictionary

            if shortest_distance != 9999:
                if start_residue[1] < end_residue[1]:  # this to order the dict so that chain goes alphabetically
                    sasds[start_residue, end_residue, shortest_distance] = all_distances[shortest_distance]
                elif end_residue[1] < start_residue[1]:
                    sasds[end_residue, start_residue, shortest_distance] = all_distances[shortest_distance]
                # if both on the same chain, then ordered to go numerically
                elif start_residue[0] < end_residue[0]:
                    sasds[start_residue, end_residue, shortest_distance] = all_distances[shortest_distance]
                else:
                    sasds[end_residue, start_residue, shortest_distance] = all_distances[shortest_distance]

    return sasds

def calculate_SASDs_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return calculate_SASDs(*a_b)

def calculate_specific_SASD_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return calculate_specific_SASD(*a_b)

def parallel_BFS(aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA, crosslink_pairs,
                 max_dist, vox, ncpus, xl_list):

    """

    Parallelised Breadth First Search of grid.

    Returns dictionary containing all solvent accessible surface distances
    {start res, end res, length in angstroms : voxel path of sasd}

    Arguments:

       *start_residue*
           key of aa1_voxels. aa1_voxels[start_residue] = all the starting voxels for that
           residue
       *aa1_voxels*
           dictionary containing starting voxels {start_residue : starting voxels}
       *aa2_voxels*
           dictionary containing ending voxels {end_residue : ending voxels}
       *dens_map*
           grid with solvent accessible surface (masked array)
       *aa1_CA*
           dictionary containing voxel of C-alpha
       *aa2_CA*
           dictionary containing voxel of C-alpha
       *crosslink_pairs*
           list of pairs of crosslinks (empty if not calculating specific crosslinks)
       *max_dist*
           maximum distance BFS will search until
       *vox*
           number of angstoms per voxel
       *ncpus*
           number of allocated cpus

    """

    freeze_support()
    final_XL = {}

    if xl_list:
        if ncpus > 1:
            pool = Pool(ncpus)
            xl_dictionaries = pool.map(calculate_specific_SASD_star,
                                       zip(crosslink_pairs,
                                       itertools.repeat(aa1_voxels),
                                       itertools.repeat(aa2_voxels),
                                       itertools.repeat(dens_map),
                                       itertools.repeat(aa1_CA),
                                       itertools.repeat(aa2_CA),
                                       itertools.repeat(max_dist),
                                       itertools.repeat(vox)))

            for c in xl_dictionaries:
                final_XL.update(c)

        else:
            # alternative call to allow single cpu running on Windows machines
            for single_crosslink in crosslink_pairs:
                xl_dictionaries = calculate_specific_SASD(single_crosslink, aa1_voxels,
                                                          aa2_voxels, dens_map, aa1_CA, aa2_CA,
                                                          max_dist, vox)
                final_XL.update(xl_dictionaries)

    else:
        if ncpus > 1:

            pool = Pool(ncpus)
            xl_dictionaries = pool.map(calculate_SASDs_star,
                                       zip(aa1_voxels,
                                       itertools.repeat(aa1_voxels),
                                       itertools.repeat(aa2_voxels),
                                       itertools.repeat(dens_map),
                                       itertools.repeat(aa1_CA),
                                       itertools.repeat(aa2_CA),
                                       itertools.repeat(max_dist),
                                       itertools.repeat(vox)))

            for c in xl_dictionaries:
                final_XL.update(c)

        else:
            # alternative call to allow single cpu running on Windows machines
            for start_residue in aa1_voxels:
                xl_dictionaries = calculate_SASDs(start_residue, aa1_voxels, aa2_voxels,
                                                  dens_map, aa1_CA, aa2_CA, max_dist, vox)
                final_XL.update(xl_dictionaries)

    return final_XL

def calculate_distance(cords):
    ''' Calculates the distance of points in 3d, input e.g. [[x1,y1,z1],[x2,y2,z3]] '''
    return math.sqrt(((cords[0][0]-cords[1][0])**2)+((cords[0][1]-cords[1][1])**2)+((cords[0][2]-cords[1][2])**2))

def get_euclidean_distances(sasds, pdb, aa1, aa2):

    residues = {}
    euc_dists = {}
    with open (pdb) as inf:
        for line in inf:
            if line.startswith('ATOM') and (line[12:16].strip() == 'CA'):
                if line[21:22].strip() == "":
                    chain = " "
                else:
                    chain = line[21:22].strip()
                residues[line[22:26].strip(),chain] = [float(line[30:38].strip()),
                float(line[38:46].strip()),
                float(line[46:54].strip())]

    for k,v in residues.items():
        for k1,v1 in residues.items():
            if k1 != k:

                euc_dists[int(k[0]),k[1], int(k1[0]),k1[1]] = calculate_distance([v,v1])

    sasds_and_eucs = {}

    for s in sasds:
        if (s[0][0],s[0][1],s[1][0],s[1][1]) in euc_dists:
            sasds_and_eucs[s[0],s[1],s[2],euc_dists[(s[0][0],s[0][1],s[1][0],s[1][1])]] = sasds[s]

    return sasds_and_eucs
