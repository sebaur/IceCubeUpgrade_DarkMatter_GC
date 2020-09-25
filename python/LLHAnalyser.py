import numpy as np
import scipy.special as sps
from iminuit import Minuit

class Profile_Analyser:

    def __init__(self):
        
        self.LLHtype = None

        self.ready = False
        self.signalPDF = None
        self.backgroundPDF = None

        self.signalPDF_uncert2 = None
        self.backgroundPDF_uncert2 = None

        self.nTotalEvents = 0.
        self.nSignalEvents = 0.
        
        self.nbins = 0
        
        self.observation = None

        self.computedBestFit = False
        self.bestFit = None
        self.TS = None
        
        self.moreOutput = False

    def setLLHtype(self,type):
        availableTypes = ['Poisson', 'Effective']
        if type not in availableTypes:
            raise ValueError('LLH type not implemented yet. Choose amongst: '+str(availableTypes))
        else:
            self.LLHtype = type
            
    def saveMoreOutput(self):
        self.moreOutput = True

    def loadBackgroundPDF(self,pdf):
        self.backgroundPDF = pdf.flatten()/np.sum(pdf)
        self.nTotalEvents = np.sum(pdf)
        print('total expected events:', self.nTotalEvents)
        
        self.nbins = len(pdf.flatten())

    def loadSignalPDF(self,pdf):
        self.signalPDF = pdf.flatten()/np.sum(pdf)
        self.nSignalEvents = np.sum(pdf)
        if self.nbins == len(pdf.flatten()):
            self.ready = True
        else:
            raise ValueError('Shape of signal pdf does not match the background pdf! Did you initialize the background pdf first?')
    
    def loadUncertaintyPDFs(self,bkg_pdf,sig_pdf):
        self.backgroundPDF_uncert2 = bkg_pdf.flatten()/self.nTotalEvents
        self.signalPDF_uncert2 = sig_pdf.flatten()/self.nSignalEvents

        if self.nbins != len(bkg_pdf.flatten()):
            raise ValueError('Shape of background uncertainty pdf does not match the background pdf!')
        if self.nbins != len(sig_pdf.flatten()):
            raise ValueError('Shape of signal uncertainty pdf does not match the background pdf!')
            
    def sampleObservation(self,xi):
        if not self.ready:
            raise ValueError('Not all pdfs are correctly loaded!')

        observationPDF = self.nTotalEvents* ((1-xi)*self.backgroundPDF + xi*self.signalPDF)

        self.observation=np.zeros(np.shape(self.backgroundPDF))
        for i in range(len(self.observation)):
            self.observation[i]=np.random.poisson(observationPDF[i])
        
        self.computedBestFit = False
        
    def evaluateLLH(self, xi):
        
        if self.LLHtype == 'Poisson':

            modelPDF = self.nTotalEvents*((1-xi)*self.backgroundPDF + xi*self.signalPDF)
            
            if np.isnan(modelPDF).any():
                print('nan in model array with xi=',xi, self.computedBestFit)

            bins_to_use = (modelPDF>0.)
            values = self.observation[bins_to_use]*np.log(modelPDF[bins_to_use])-modelPDF[bins_to_use]
            
            
        elif self.LLHtype == 'Effective':
        
            modelPDF = self.nTotalEvents*((1-xi)*self.backgroundPDF + xi*self.signalPDF)
            modelPDF_uncert2 = self.nTotalEvents*((1-xi)*(1-xi)*self.backgroundPDF_uncert2 + xi*xi*self.signalPDF_uncert2)

            if np.isnan(modelPDF).any():
                print('nan in model array with xi=',xi, self.computedBestFit)
                
            bins_to_use = (modelPDF>0.)&(modelPDF_uncert2>0.) 

            alpha = modelPDF[bins_to_use]**2/modelPDF_uncert2[bins_to_use] +1.
            beta  = modelPDF[bins_to_use]/modelPDF_uncert2[bins_to_use]

            values = [
              alpha*np.log(beta),
              sps.loggamma(self.observation[bins_to_use]+alpha).real,
              -(self.observation[bins_to_use]+alpha)*np.log1p(beta),
              -sps.loggamma(alpha).real,
            ]
            
        else:
            raise ValueError('No valid LLH type defined!')
        
        return -np.sum(values)

    
    def ComputeBestFit(self):
        
        lower_bound = 0.
        LLHmin_DM=Minuit(self.evaluateLLH,
             xi=0.1, error_xi=.01,
             limit_xi=(lower_bound,2.),
             errordef=.5,print_level=0)  
        LLHmin_DM.migrad()
        
        self.bestFit = {}
        self.bestFit['xi']=LLHmin_DM.fitarg['xi']
        self.bestFit['LLH']=self.evaluateLLH(self.bestFit['xi'])
                
        self.computedBestFit = True
        
        
    def ComputeTestStatistics(self):
        
        if not self.computedBestFit:
            self.ComputeBestFit()

        self.TS = 0.
        self.bestFit['LLH_ref'] = self.evaluateLLH(0.)
        
        if self.bestFit['xi'] > 0.:
            self.TS = 2*(self.bestFit['LLH_ref']-self.bestFit['LLH'])
        

    def CalculateUpperLimit(self,conf_level):

        nIterations = 0
        eps_TS=0.005
        eps_param=0.0005

        deltaTS = 2.71
        if conf_level==90:
            deltaTS = 1.64
        elif conf_level==95:
            deltaTS = 2.71
            
        param_low=self.bestFit['xi']
        param_up=self.bestFit['xi']
        param_mean=self.bestFit['xi']
        
        dTS=0
        cc=1
        while((dTS<deltaTS) and (nIterations<100)):
            nIterations += 1 

            param_up=param_up+3.*np.abs(param_up)

            if param_up <0.:
                TS_fix = 0.
            else:
                TS_fix = 2*(self.bestFit['LLH_ref']-self.evaluateLLH(param_up))
                            
            dTS = self.TS - TS_fix

        nIterations = 0
        param_low=param_up/4.
        while((cc>0.)  and (nIterations<100)):
            
            nIterations += 1

            param_mean=(param_low+param_up)/2.
            if param_mean <0.:
                TS_fix = 0.
            else:
                TS_fix = 2*(self.bestFit['LLH_ref']-self.evaluateLLH(param_mean))
                
            dTS = self.TS - TS_fix
            
            if(dTS<deltaTS):

                param_low=param_mean
                delta_param=(param_up-param_low)/(param_up)
                
                if((dTS>deltaTS-eps_TS) and (delta_param < eps_param)):
                    cc = 0
                    
            if(dTS>deltaTS):
                param_up=param_mean
                delta_param=(param_up-param_low)/(param_up)
                
                if((dTS<deltaTS+eps_TS) and (delta_param < eps_param)):
                    cc=0
                    
        return param_up
       
    
    def CalculateSensitivity(self,nTrials, conf_level):

        if self.LLHtype == None:
            raise ValueError('LLH type not defined yet!')

        TS = []
        upperlimits = []
        if self.moreOutput:
            fits = []
        
        for i in range(nTrials):
            self.sampleObservation(0.)
            self.ComputeTestStatistics()
            TS.append(self.TS)

            ul = self.CalculateUpperLimit(conf_level)
            if np.isnan(ul):
                print("Warning: NaN upper limit at trial {i}.\nRepeating trial.".format(i=i))
                i-=1
                continue
            upperlimits.append(ul)
            
            if self.moreOutput:
                fits.append(self.bestFit)
            
        p_median = np.percentile(upperlimits, 50)
        p_95_low = np.percentile(upperlimits, 2.5)
        p_95_high = np.percentile(upperlimits, 97.5)
        p_68_low = np.percentile(upperlimits, 16.)
        p_68_high = np.percentile(upperlimits, 84.)

        dic_brazilian = {}
        dic_brazilian['TS_dist'] = TS
        dic_brazilian['error_68_low'] = p_68_low
        dic_brazilian['error_68_high'] = p_68_high
        dic_brazilian['error_95_low'] = p_95_low
        dic_brazilian['error_95_high'] = p_95_high   
        dic_brazilian['median'] = p_median
        if self.moreOutput:
            dic_brazilian['upperlimits'] = upperlimits
            dic_brazilian['bestFits'] = fits

        return dic_brazilian   
