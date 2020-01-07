import scipy.stats as sts

class Combine:
    def __init__(self,MNXL,MoDS):
        self.mnxl = MNXL
        self.mods = MoDS

        self.zmnxl, self.zmods = self.normalise_scores()
        self.zxlmo = self.score_xlmo()

    def normalise_scores(self):
        return sts.zscore(self.mnxl), sts.zscore(self.mods)

    def score_xlmo(self):
        return sts.zscore([i+j for i,j in zip(self.zmnxl,self.zmods)])
