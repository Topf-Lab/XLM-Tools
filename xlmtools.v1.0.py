import os
import sys
from XLMTools import MNXL, MoDS, XLMO
import argparse
import subprocess as sp
import Jwalk.RunJwalk as rj
from multiprocessing import cpu_count

def run_depth(pdbs,depth_source='DEPTH'):
    try:
        os.mkdir("depth_files")
    except:
        pass

    try:
        for pdb in pdbs:
            cmd = [depth_source,"-i",pdb,"-o","depth_files/%s" % pdb.split(".pdb")[0]]
            sp.call(cmd)
            #test to see if file created, if not will give error
            with open("depth_files/%s-residue.depth" % pdb.split(".pdb")[0]) as f:
                f.close()
    except:
        print("Residue depth calculation failed. Please check DEPTH installation or specify depth source (see help, -h flag).")
        print("Please note residue depth calculation may not work on all OS. Please see http://cospi.iiserpune.ac.in/depth/htdocs/download.html for details.")
        sys.exit(0)

jwalk_files = []
depth_files = []
sep = "\t"

parser = argparse.ArgumentParser(prog="XLM-Tools",description="XLM Tools tools: A tool to score model protein structures according to crosslink and monolink data.")
xlm_specific = parser.add_argument_group(title="XLM Arguments")
xlm_specific.add_argument('-xl_list',nargs=1,
                    help="input list of crosslinks and monolinks")
xlm_specific.add_argument('-jwalk_files',nargs="*",
                    help="jwalk files for MNXL scoring")
xlm_specific.add_argument('-depth_files',nargs='+',
                    help="depth files for MoDS scoring")
xlm_specific.add_argument('-outfile_name',nargs=1,
                    help='specify output file name')
xlm_specific.add_argument('-sep', nargs=1,
                    help="separator in output file, default is tab")
xlm_specific.add_argument('-pdb',nargs='+',
                    help="specify pdb files for Jwalk/Depth run")
#commands for running jwalk
jwalk_args = parser.add_argument_group(title="Jwalk Arguments")
jwalk_args.add_argument('-jwalk',action="store_true",
                    help='flag to use if starting from .pdb files and running Jwalk')
jwalk_args.add_argument('-vox', nargs=1,
					help='specify voxel size of grid')
jwalk_args.add_argument('-surface', action="store_true",
					help='use higher accuracy method to calculate solvent accessibility - requires Freesasa installation')
jwalk_args.add_argument('-ncpus', nargs=1,
					help='specify number of cpus to use')
#commands for running depth
depth_args = parser.add_argument_group(title="Depth Arguments")
depth_args.add_argument('-depth',action='store_true',
                        help="flag to use if starting from .pdb files and running Jwalk")
depth_args.add_argument('-depth_source',nargs=1,
                        help='specify depth source')


args = parser.parse_args()

if args.xl_list:
    xl_list = args.xl_list[0]
else:
    raise Exception('Please provide Xl-List')

if args.pdb:
    pdb_list = args.pdb
else:
    pdb_list = [i for i in os.listdir('.') if i.endswith('.pdb')]

if args.jwalk:
    if args.vox:
        vox = args.vox[0]
    else:
        vox = 1

    if args.surface:
        surface = True
    else:
        surface = False

    if args.ncpus:
        ncpus = args.ncpus[0]
    else:
        ncpus = cpu_count()

    rj.runJwalk(100,vox,surface,False,False,xl_list,"LYS","LYS",ncpus,pdb_list)
    jwalk_files = ["Jwalk_results/%s" % i for i in os.listdir('Jwalk_results/') if i.endswith('_crosslink_list.txt')]

if args.jwalk_files:
    jwalk_files = args.jwalk_files
else:
    if not jwalk_files:
        jwalk_files = [i for i in os.listdir('.') if i.endswith('_crosslink_list.txt')]

if args.depth:
    if args.depth_source:
        depth_source = args.depth_source[0]
    else:
        depth_source = 'DEPTH'

    run_depth(pdb_list)
    depth_files = ["depth_files/%s" % i for i in os.listdir('depth_files') if i.endswith("-residue.depth")]

if args.depth_files:
    depth_files = args.depth_files
else:
    if not depth_files:
        depth_files = [i for i in os.listdir('.') if i.endswith('-residue.depth')]

if args.outfile_name:
    outfile_name = args.outfile_name
else:
    outfile_name = "xlm_scores.txt"

if args.sep:
    sep = args.sep[0]

class Reference:
    def __init__(self, xl_list, jwalk_files=None, depth_files=None):
        #Initialise variables
        self.ml = []
        self.intra = []
        self.inter = []

        self.chains = []


        #if neither provided need to run jwalk/depth
        self.jwalk_files = jwalk_files
        self.depth_files = depth_files

        self.results = {i.split('_crosslink')[0].split('/')[-1]:{} for i in jwalk_files}

        #opens xl_list, parses for xl/ml
        with open(xl_list) as f:
            g = f.read().splitlines()
            for line in g:
                col = line.split('|')
                if len(col) == 4:
                    aa1 = int(col[0])
                    c1 = col[1]
                    aa2 = int(col[2])
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

                    if c1 not in self.chains:
                        self.chains.append(c1)
                    if c2 not in self.chains:
                        self.chains.append(c2)

                else:
                    aa = int(col[0])
                    c = col[1]
                    self.ml.append([aa,c])
            f.close()

        if len(self.chains) > 1:
            self.complex = 1
        else:
            self.complex = 0

        """
        Control flow dictates what calculations are performed:

        if only XL -> MNXL
        if only ML -> MoDS
        if both -> MNXL, MoDS and XLMO
        """

        if self.ml and not self.intra:
            print('Calculating MoDS on all files...')
            self.out_option = 'MoDS'
            if self.depth_files:
                for d in depth_files:
                    self.mods_only(d)
            else:
                raise Exception("Depth files missing, please specify path.")
        elif self.intra and not self.ml:
            print('Calculating MNXL on all files...')
            self.out_option = 'MNXL'
            if self.jwalk_files:
                for j in jwalk_files:
                    self.mnxl_only(j)
            else:
                raise Exception("Jwalk files missing, please specify path.")
        elif self.intra and self.ml:
            print('Calculating MNXL, MoDS and XLMO on all files...')
            self.out_option = 'XLMO'
            if self.depth_files and self.jwalk_files:
                for j in self.jwalk_files:
                    jwalk = j
                    if jwalk.startswith("Jwalk_results"):
                        depth = "depth_files/%s-residue.depth" % j.split('_crosslink')[0].split('/')[-1]
                    else:
                        depth = "%s-residue.depth" % j.split('_crosslink')[0]
                    self.mnxl_only(jwalk)
                    print(self.results)
                    self.mods_only(depth)
            elif not self.depth_files:
                raise Exception("Depth files missing, please specify path.")
            else:
                raise Exception("Jwalk files missing, please specify path.")

        if self.out_option == 'XLMO':
            self.score_xlmo()

    def mnxl_only(self, file):
        m = file

        file_pass = False
        with open(m) as check_file:
            if next(check_file)[0:5] == "Index":
                file_pass = True

        if file_pass:
            # load the model data in and score
            # print(m)
            model = MNXL.Model(m)
            model.load_crosslinks()
            # print(model.intra)
            model.score_intra(self)
            model.score_inter(self)
            model.penalise_missing_intra(self)
            model.penalise_missing_inter(self)
            model.generate_totals()

        if self.complex == 0:
            self.results[m.split('_crosslink_list.txt')[0].split('/')[-1]].update({'MNXL':model.MNXL, 'Matched':model.number_of_matched, 'NoV':model.number_of_violations, 'NoNA':model.number_of_non_access})
        else:
            self.results[m.split('_crosslink_list.txt')[0].split('/')[-1]].update({'cMNXL':model.cMNXL, 'Matched':model.number_of_matched, 'NoV':model.number_of_violations, 'NoNA':model.number_of_non_access})

    def mods_only(self, file):
        print(file)
        model = MoDS.Depth(file)
        model.load_monolinks()
        model.score_mono(self)
        self.results[file.split('-residue.depth')[0].split('/')[-1]].update({'MoDS':model.mods})

    def score_xlmo(self):

        if 'cMNXL' in self.results[list(self.results.keys())[0]]:
            mnxl_key = 'cMNXL'
        else:
            mnxl_key = 'MNXL'

        names = []
        mnxl = []
        mods = []

        for k,v in self.results.items():
            names.append(k)
            mnxl.append(v[mnxl_key])
            mods.append(v['MoDS'])

        xlmo = XLMO.Combine(mnxl,mods)

        # self.zmnxl, self.zmods = xlmo.normalise_scores()
        xlmo_scores = xlmo.score_xlmo()

        for name,score,zmn,zmo in zip(names,xlmo_scores, xlmo.zmnxl, xlmo.zmods):
            self.results[name].update({'XLMO':score,'zMNXL':zmn,'zMoDS':zmo})


    def output_scores(self,name="xlm_scores.txt",sep="\t"):
        with open(name,"w+") as f:
            titles = [i for i in self.results[list(self.results.keys())[0]]]
            f.write("Model%s%s\n" % (str(sep),str(sep).join(titles)))

            for k,v in self.results.items():
                # print(k)
                values = []
                for t in titles:
                    values.append(v[t])
                f.write("%s%s%s\n" % (k.split('/')[-1],str(sep),str(sep).join([str(i) for i in values])))
            f.close()

        print("\nScores written to %s. Please cite Sinnott et al., Combining Information from Crosslinks and Monolinks in the Modelling of Protein Structures, Structure (2020), https://doi.org/10.1016/j.str.2020.05.012" % name)






######################################################

if __name__ == "__main__":
    ref = Reference(xl_list,jwalk_files,depth_files)
    ref.output_scores(name=outfile_name,sep=sep)
