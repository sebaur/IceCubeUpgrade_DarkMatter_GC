import pandas as pd
import math

from oscillations import infOsz
import scipy
from scipy.interpolate import griddata, PchipInterpolator

import numpy as np

pdg2flav = {}
pdg2flav[12] = 'nuE'
pdg2flav[14] = 'nuMu'
pdg2flav[16] = 'nuTau'
pdg2flav[-12] = 'nuEBar'
pdg2flav[-14] = 'nuMuBar'
pdg2flav[-16] = 'nuTauBar'

class PPPC_spectra():

    def __init__(self,path, mass, process):
        self.basePath = path
        self.mass = mass
        self.process = process
        self.columns = {'nuE': '\[Nu]e', 'nuMu': '\[Nu]\[Mu]', 'nuTau': '\[Nu]\[Tau]', 'mu': '\[Mu]', 'b': 'b', 'W': 'W', 'tau': '\[Tau]'}
        self.loaded = False
        
        self.vec_getSourceSpectrumWeight = np.vectorize(self.getSourceSpectrumWeight, excluded='self')
        self.vec_getEarthSpectrumWeight = np.vectorize(self.getEarthSpectrumWeight, excluded='self')


    def load(self):
        self.rawData = {}
        try:
            self.rawData['nuE'] = pd.read_csv(self.basePath+'/AtProduction_neutrinos_e.dat', sep=r"\s*", engine='python')
            self.rawData['nuMu'] = pd.read_csv(self.basePath+'/AtProduction_neutrinos_mu.dat', sep=r"\s*", engine='python')
            self.rawData['nuTau'] = pd.read_csv(self.basePath+'/AtProduction_neutrinos_tau.dat', sep=r"\s*", engine='python')
            self.loaded = True
        except:
            raise NameError('PPPC spectra not found at path ' + self.basePath +' !')

        self.log10x = self.rawData['nuE'].loc[(self.rawData['nuE'].mDM == self.mass), ['Log[10,x]']].values.T[0]

        self.sourceSpectra = {}
        self.sourceSpectra['nuE'] = self.rawData['nuE'].loc[(self.rawData['nuE'].mDM == self.mass), [self.columns[self.process]]].values.T[0]
        self.sourceSpectra['nuMu'] = self.rawData['nuMu'].loc[(self.rawData['nuMu'].mDM == self.mass), [self.columns[self.process]]].values.T[0]
        self.sourceSpectra['nuTau'] = self.rawData['nuTau'].loc[(self.rawData['nuTau'].mDM == self.mass), [self.columns[self.process]]].values.T[0]

        oszillation = infOsz()
        oszillatedSpectra = oszillation.oscillate(self.sourceSpectra['nuE'],self.sourceSpectra['nuE'],
                                    self.sourceSpectra['nuMu'],self.sourceSpectra['nuMu'],
                                    self.sourceSpectra['nuTau'],self.sourceSpectra['nuTau'])
        self.earthSpectra = {}
        self.earthSpectra['nuE'] = oszillatedSpectra[0] / 2.
        self.earthSpectra['nuMu'] = oszillatedSpectra[1] / 2.
        self.earthSpectra['nuTau'] = oszillatedSpectra[2] / 2.

        self.sourceSpectrum_interpol = {}
        self.earthSpectrum_interpol = {}

        for f in ['nuE','nuMu','nuTau']:
            self.earthSpectrum_interpol[f] = PchipInterpolator(
                self.log10x, self.earthSpectra[f], extrapolate=False)    
            self.sourceSpectrum_interpol[f] = PchipInterpolator(
                self.log10x, self.sourceSpectra[f], extrapolate=False)    
        
    def getSourceSpectrum(self,flavor):
        if not self.loaded:
            raise NameError('PPPC spectra not loaded!')   
        return (self.log10x, self.sourceSpectra[flavor])

    def getEarthSpectrum(self,flavor):
        if not self.loaded:
            raise NameError('PPPC spectra not loaded!')   
        return (self.log10x, self.earthSpectra[flavor])

    def getEarthSpectrumWeight(self,pdg,E):
        x = E/self.mass
        return self.earthSpectrum_interpol[pdg2flav[np.abs(pdg)]](np.log10(x))/(np.log(10)*x)/self.mass

    def getSourceSpectrumWeight(self,pdg,E):
        x = E/self.mass
        return self.sourceSpectrum_interpol[pdg2flav[np.abs(pdg)]](np.log10(x))/(np.log(10)*x)/self.mass
