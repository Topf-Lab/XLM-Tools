import sys
import os
import numpy as np
import argparse
import math
import numpy
import re
import scipy.stats
import subprocess, shlex

"""
MNXL needs as input:

1) reference data (i.e. the experimental crosslinks) formatted as shown in mnxl_reference_example.txt.

2) At least one, but ideally a list of Jwalk output files to score
OR
2) pdb files that can be passed via Jwalk to create output files.

Other than calling jwalk externally, everything else is contained in this script

"""

class Reference:
    """Object to handle the reference data"""

    def __init__(self,ref_xl):

        chains = []

        with open(ref_xl) as f:
            g = f.read().splitlines()

            self.intra = []
            self.inter = []

            for line in g:
                l = line.rstrip()
                col = l.split('|')
                if len(col) == 4:
                    aa1 = col[0]
                    c1 = col[1]
                    aa2 = col[2]
                    c2 = col[3]
                    if c1 == c2:
                        if int(aa1) < int(aa2):
                            low = aa1; low_chain = c1; hi = aa2; hi_chain = c2
                        else:
                            low = aa2; low_chain = c2; hi = aa1; hi_chain = c1

                        self.intra.append([low,low_chain,hi,hi_chain])

                    else:
                        if int(aa1) < int(aa2):
                            low = aa1; low_chain = c1; hi = aa2; hi_chain = c2
                        else:
                            low = aa2; low_chain = c2; hi = aa1; hi_chain = c1

                        self.inter.append([low,low_chain,hi,hi_chain])

                    if c1 not in chains:
                        chains.append(c1)
                    if c2 not in chains:
                        chains.append(c2)



        if len(chains) > 1:
            self.complex = 1
        else:
            self.complex = 0

class Model():
    """Object to handle each model's data"""

    def __init__(self,file,weight=0.3,cut=36):

        self.file = file
        # file_search = re.search('native',file,re.IGNORECASE)
        self.inter_cut = cut
        self.intra_weight = weight
        self.fnum = file

        self.P_inter = 0
        self.P_intra = 0
        self.N_inter = 0
        self.N_intra = 0
        self.M_inter = 0
        self.M_intra = 0

        self.number_of_violations = 0
        self.number_of_matched = 0
        self.number_of_non_access = 0

    def load_crosslinks(self):
        '''self.intra, self.inter'''
        with open(self.file) as f:
            next(f)
            self.intra = {};self.inter = {}
            self.euc_intra = {}; self.euc_inter = {}
            for xl in [i for i in f if i != '\n']:
                col = xl.split()
                col2 = col[2]
                col3 = col[3]

                SASD = float(col[4])
                aa1 = col2.split('-')
                aa2 = col3.split('-')
                if aa1[2] == aa2[2]:
                    if int(aa1[1]) < int(aa2[1]):
                        low = aa1[1]; low_chain = aa1[2]; hi = aa2[1]; hi_chain = aa2[2]
                    else:
                        low = aa2[1]; low_chain = aa2[2]; hi = aa1[1]; hi_chain = aa1[2]
                    self.intra[int(low),low_chain,int(hi),hi_chain] = SASD
                else:
                    if int(aa1[1]) < int(aa2[1]):
                        low = aa1[1]; low_chain = aa1[2]; hi = aa2[1]; hi_chain = aa2[2]
                    else:
                        low = aa2[1]; low_chain = aa2[2]; hi = aa1[1]; hi_chain = aa1[2]

                    self.inter[int(low),low_chain,int(hi),hi_chain] = SASD

    def score_intra(self,Reference):
        '''self.intra_score, self.intra_count, self.intra_viol, self.intra_sum_viol'''
        N=scipy.stats.norm(18.62089,5.995381)

        for xl1,c1,xl2,c2 in self.intra:
            if [xl1,c1,xl2,c2] in Reference.intra:
                SASD = self.intra[xl1,c1,xl2,c2]
                if SASD <= 33:
                    self.P_intra += N.pdf(SASD)
                    self.number_of_matched += 1
                else:
                    self.N_intra += -0.1
                    self.number_of_violations += 1

    def score_inter(self,Reference):
        '''self.inter_score, self.inter_count, self.inter_viol, self.inter_sum_viol'''
        N=scipy.stats.norm(21.91883,4.871774)

        for xl1,c1,xl2,c2 in self.inter:
            if [xl1,c1,xl2,c2] in Reference.inter:
                SASD = self.inter[xl1,c1,xl2,c2]
                if SASD <= 36:
                    self.P_inter += N.pdf(SASD)
                    self.number_of_matched += 1
                else:
                    self.N_inter += -0.1
                    self.number_of_violations += 1

    def penalise_missing_intra(self,Reference):

        for ref_xl1,c1 ,ref_xl2,c2 in Reference.intra:
            if (ref_xl1,c1,ref_xl2,c2) not in self.intra:
                self.M_intra += -0.1
                self.number_of_non_access += 1

    def penalise_missing_inter(self,Reference):

        for ref_xl1,c1, ref_xl2,c2 in Reference.inter:
            if (ref_xl1,c1,ref_xl2,c2) not in self.inter:
                self.M_inter += -0.1
                self.number_of_non_access += 1

    def generate_totals(self):
        self.inter_NPM = self.N_inter+self.P_inter+self.M_inter
        self.cMNXL = self.inter_NPM + self.M_intra*0.3
        self.MNXL = self.N_intra + self.M_intra + self.P_intra

def output_scores(data,titles,outfile):
    with open(outfile, "w") as f:

        f.write('\t'.join('{}'.format(col) for col in titles))
        f.write('\n')
        for line in data:
            if line != '':
                f.write('\t'.join('{}'.format(i) for i in line))
                f.write('\n')

def calculate_cMNXL(referee,model_crosslinks):

    """
    calculates MNXL or cMNXL using the refernce data (referee) and jwalk_output_files (model_crosslinks)
    outputs score file in .txt format
    """

    ref = Reference(referee)
    data = []
    for m in model_crosslinks:
        # check Jwalk output files are formatted correctly
        file_pass = False
        with open(m) as check_file:
            if next(check_file)[0:5] == "Index":
                file_pass = True

        if file_pass:
            # load the model data in and score
            # print(m)
            model = Model(m)
            model.load_crosslinks()
            model.score_intra(ref)
            model.score_inter(ref)
            model.penalise_missing_intra(ref)
            model.penalise_missing_inter(ref)
            model.generate_totals()

            # cMNXL calculates if the model is a complex or monomer and scores appropriately
            if ref.complex == 1:
                data.append([model.fnum,
                    model.cMNXL,
                    model.number_of_matched,
                    model.number_of_violations,
                    model.number_of_non_access])
            else:
                data.append([model.fnum,
                    model.MNXL,
                    model.number_of_matched,
                    model.number_of_violations,
                    model.number_of_non_access])

    if ref.complex == 1:
        titles =  ["model",
                   "cMNXL (protein complex)",
                   "Matched",
                   "Violating",
                   "Non-accessible"]
        out_file = "cMNXL_scores.txt"
        output_scores(data,titles, out_file)
    else:
        titles =  ["model",
                   "MNXL (protein monomer)",
                   "Matched",
                   "Violating",
                   "Non-accessible"]
        out_file = "MNXL_scores.txt"
        output_scores(data,titles,out_file)

def run_jwalk(pdb_list, referee):
    """
    calls jwalk externally
    """
    for pdb in pdb_list:
        arg_string = "jwalk -max_dist 500 -i %s -xl_list %s" % (pdb, referee)
        jwalk_run_options = shlex.split(arg_string)
        result = subprocess.Popen(jwalk_run_options)
        result.communicate()


# if __name__ == "__main__":
#
#     parser = argparse.ArgumentParser(description='MNXL: Model validation using crosslinking restraints ')
#     parser.add_argument('-data', nargs=1,
#                         help='specify experimental data: -data <input_data.txt>')
#     parser.add_argument('-mod_xl', nargs="+",
#                         help='specify simulated model data: -mod_xl <model_data.txt>')
#     parser.add_argument('-jwalk', action="store_true",
#                         help='flag to use if starting from .pdb files and running Jwalk')
#     parser.add_argument('-pdb', nargs="+",
#                         help='specify input pdbs: -mod_xl <model_data.pdb>')
#
#     args = parser.parse_args()
#
#     if args.data:
#         referee = args.data[0]
#     else:
#         print "Please specify the experimental data file, use -h flag for help"
#         sys.exit()
#
#     if args.jwalk:
#
#         if args.pdb:
#             pdb_list = args.pdb
#         else:
#             pdb_list = []
#             for m in [i for i in os.listdir("./") if i.endswith(".pdb")]:
#                 pdb_list.append(m)
#
#             if len(pdb_list) < 1:
#                 print "Please specify .pdb files or place in current directory. use -h flag for help"
#                 sys.exit()
#
#         run_jwalk(pdb_list, referee)
#
#         model_crosslinks = []
#         for model_xl in [i for i in os.listdir("./Jwalk_results/") if i.endswith("list.txt")]:
#             model_crosslinks.append("./Jwalk_results/"+model_xl)
#
#     elif args.mod_xl:
#         model_crosslinks = args.mod_xl
#     else:
#         model_crosslinks = []
#         for m in [i for i in os.listdir("./") if i.endswith(".txt")]:
#             model_crosslinks.append(m)
#
#         if len(model_crosslinks) < 1:
#             print "Please specify Jwalk output files. use -h flag for help"
#             sys.exit()
#
#     calculate_cMNXL(referee, model_crosslinks)
