import pandas as pd
import math
import numpy as np

class Honda_Flux():

    def __init__(self,path):
        self.basePath = path
        self.columns = {'E': 'Enu(GeV)', 'nuE': 'NuE', 'nuEbar': 'NuEbar', 'nuMu': 'NuMu','nuMubar': 'NuMubar'}
        self.loaded = False


    def load(self):
        self.rawData = {}
        try:
            for zenithBin in range(20):
                self.rawData[zenithBin] = pd.read_csv(self.basePath+'/Honda_SP.dat', sep=r"\s*", engine='python',skiprows=(zenithBin*102+(zenithBin+1)*1),nrows=101)
            self.loaded = True
        except:
            raise NameError('Honda tables not found at path ' + self.basePath +' !')

    def getFlux(self, flavor, zenithBin):
        if not self.loaded:
            raise NameError('Honda tables not loaded!')   
        return (self.rawData[zenithBin].loc[: , [self.columns["E"]]].values.T[0],
                self.rawData[zenithBin].loc[: , [self.columns[flavor]]].values.T[0])
