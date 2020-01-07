class Reference:
    """A class to store the monolink information from an XLMS input file.
    XL info will be dealtwith in MNXL.py"""
    def __init__(self, ref_ml):
        self.mono = []
        with open(ref_ml) as f:
            g = f.read().splitlines()
            for line in [i for i in g if len(i.split('|')) == 2]:
                col = line.split('|')
                aa = col[0]
                chain = col[1]
                self.mono.append((int(aa),chain))
            f.close()

class Depth:
    """Object to handle depth information"""
    def __init__(self,file):
        self.file = file
        self.mods = 0

    def load_monolinks(self):
        self.mono = {}
        with open(self.file) as f:
            g = f.read().splitlines()
            for line in [i for i in  g if not i.startswith('#')]:
                mono_id = int(line.split()[0].split(':')[1])
                mono_chain = line.split()[0].split(':')[0]
                self.mono[mono_id,mono_chain] = float(line.split()[2])
            f.close()

    def score_mono(self,Reference):
        for ml,c in self.mono:
            if [ml,c] in Reference.ml:
                dp = self.mono[(ml,c)]
                self.mods-=dp
