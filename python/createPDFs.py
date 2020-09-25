#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /data/user/sbaur/metaprojects/icerec/build/

import time
from time import gmtime, strftime
import numpy as np
import sys, os, copy
import pickle
from optparse import OptionParser
import random
from icecube import astro
from astropy.time import Time as apTime
import pandas as pd
from scipy import interpolate

#######################
# get and define parameters
#######################

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-t", "--type",default="background",
                  dest="TYPE", help="Define type of PDF")
parser.add_option("-c", "--channel",default="nu",
                  dest="CHANNEL", help="Annihilation channel")
parser.add_option("-p", "--profile",default="NFW",
                  dest="PROFILE", help="DM Profile")
parser.add_option("-m", "--mass", default='0',
                  dest="MASS", help="DM mass")
parser.add_option("-b", "--binning", default='RA-DEC',
                  dest="BINNING", help="Bins for PDF")
parser.add_option("-o", "--oversampling", default='100',
                  dest="OVERSAMPLING", help="Oversampling")


(options,args) = parser.parse_args()
        
mode = options.TYPE #annihilation, decay, background
channel = options.CHANNEL
profile = options.PROFILE
binning = options.BINNING
nOversampling = int(options.OVERSAMPLING)

if binning not in ['RA-DEC','Psi-E']:
    print 'Binning not defined!'
    sys.exit(1)
    
if mode != 'background':
    mass = float(options.MASS)

source = 'source /data/user/sbaur/projects/PrepareForGit/Upgrade_DarkMatter/env.sh'
dump = '/usr/bin/python -c "import os,pickle;print pickle.dumps(os.environ)"'
penv = os.popen('%s && %s' %(source,dump))
env = pickle.loads(penv.read())
os.environ = env

base_path = os.environ['ANALYSIS_BASE_PATH']
out_path = base_path+'/PDFs/'+mode+'/'

flavors = ['nuE','nuMu','nuTau']
GC_Pos = (266.*np.pi/180. , -29.*np.pi/180)


sys.path.append(base_path+'/resources')
print sys.path

def oversample(tmp_weight, tmp_nu_type,
               tmp_energy_true, tmp_energy_reco,
               tmp_zenith_true, tmp_zenith_reco,
               tmp_azimuth_true, tmp_azimuth_reco,
               nOversampling):
    
    stime = time.mktime(time.strptime("1/1/2022/00/00/00", '%m/%d/%Y/%H/%M/%S'))
    etime = time.mktime(time.strptime("1/1/2023/00/00/00", '%m/%d/%Y/%H/%M/%S'))
    
    
    oversampled_psi_true = np.zeros((nOversampling,len(tmp_weight)))
    oversampled_psi_reco = np.zeros((nOversampling,len(tmp_weight)))
    
    oversampled_RA_reco = np.zeros((nOversampling,len(tmp_weight)))
    oversampled_Dec_reco = np.zeros((nOversampling,len(tmp_weight)))


    for i in range(nOversampling):

        eventTime = stime + random.random() * (etime - stime)
        eventTime = apTime(eventTime,format='unix').mjd

        RA_true, DEC_true = astro.dir_to_equa(tmp_zenith_true,tmp_azimuth_true,eventTime)
        RA_reco, DEC_reco = astro.dir_to_equa(tmp_zenith_reco,tmp_azimuth_reco,eventTime)

        oversampled_psi_true[i] = astro.angular_distance(RA_true,DEC_true,GC_Pos[0],GC_Pos[1])
        oversampled_psi_reco[i] = astro.angular_distance(RA_reco,DEC_reco,GC_Pos[0],GC_Pos[1])

        oversampled_RA_reco[i] = RA_reco
        oversampled_Dec_reco[i] = DEC_reco
        
        
    # append n versions
    oversampled_weight = np.tile(tmp_weight/nOversampling,nOversampling)
    oversampled_nu_type = np.tile(tmp_nu_type,nOversampling)
    oversampled_energy_true = np.tile(tmp_energy_true,nOversampling)
    oversampled_energy_reco = np.tile(tmp_energy_reco,nOversampling)

    return oversampled_weight, oversampled_nu_type, oversampled_energy_true, oversampled_energy_reco, oversampled_psi_true.flatten(), oversampled_psi_reco.flatten(), oversampled_RA_reco.flatten() , oversampled_Dec_reco.flatten()


out_file = ''

if mode == 'signal':
    out_file = out_path+'Signal_PDF_'+profile+'_'+str(int(mass))+'GeV_'+str(channel)+'_oversampling'+str(nOversampling)+'_'+binning+'.pkl'
    
if mode == 'background':
    out_file = out_path+'Background_PDF_Honda2015_oversampling'+str(nOversampling)+'_'+binning+'.pkl'
    
if os.path.isfile(out_file):
    print " PDF file exists already -- won't do anything!"
    sys.exit(0)
    
    
#######################
# load Upgrade dataset
#######################

print 'load Upgrade MC'

input_file = os.environ['DATA_RELEASE_DIR']+'/events/neutrino_mc.csv'
input_data = pd.read_csv(input_file)

print '...done'

if mode == 'signal':
    #######################
    # load PPPC spectra
    #######################

    print 'load PPPC spectra'
    from PPPC4 import PPPC_spectra

    if mass < 5.:
        spectra = PPPC_spectra(base_path+'/resources/PPPC_GC/',5,channel)
        spectra.load()

    elif channel=='b' and mass == 5.:
        print 'this is a special case!'
        spectra = PPPC_spectra(base_path+'/resources/PPPC_GC/',6,channel)
        spectra.load()

    else:
        spectra = PPPC_spectra(base_path+'/resources/PPPC_GC/',mass,channel)
        spectra.load()

    print '...done'

    #######################
    # load DM profile (J)
    #######################

    print 'load J factors'

    from  Jfactors import Jfactor

    psi_sample = np.linspace(0,np.pi,1000) 

    exponent = 2
    J = interpolate.UnivariateSpline(psi_sample,Jfactor(profile, psi_sample, exponent),k=1,s=0)

    print '...done'

    #######################
    # look at data
    #######################

    print 'process data'

    mask = input_data['true_energy']<mass
    simWeight = np.array(input_data['weight'][mask], dtype=float)

    spectrum_weight = np.array(
        spectra.vec_getEarthSpectrumWeight(
            np.array(input_data['pdg'][mask]),
            np.array(input_data['true_energy'][mask])
        )
    )
    
    oversampled_weight, oversampled_pdg, oversampled_energy_true, oversampled_energy_reco, oversampled_psi_true, oversampled_psi_reco, oversampled_RA_reco, oversampled_Dec_reco = oversample(simWeight*spectrum_weight, np.array(input_data['pdg'][mask]),input_data['true_energy'][mask], input_data['reco_energy'][mask], input_data['true_zenith'][mask], input_data['reco_zenith'][mask], input_data['true_azimuth'][mask], input_data['reco_azimuth'][mask], nOversampling)

    Jfactor_weight = np.array(J(oversampled_psi_true))

    DMflux_weight = np.array(1./(2*4*np.pi*mass**2) * Jfactor_weight,dtype=float) *1E4 # from cm^-2 to m^-2
    totalWeight = oversampled_weight * DMflux_weight

    print '...done'
    
    
if mode == 'background':
    #######################
    # load Honda background
    #######################

    from HondaFlux import Honda_Flux
    from oscillations import threeFlavorVacuumOsc

    print 'load Honda flux'

    Honda = Honda_Flux(base_path+'/resources/Honda/')
    Honda.load()

    def getHondaWeight(E,cosZenith,pdg):
        cosZenithBins = np.linspace(1, -1, 21)
        iEnergyBin = np.searchsorted(Honda.getFlux('nuE',0)[0], E)
        iCosZenithBin = np.searchsorted(-cosZenithBins, -cosZenith)

        earth_radius = 6371.
        production_height = 15. # Assuming neutrino produced 15 km above surface
        detector_depth = 1. # Assuming detector depth of 1 km
        baseline = -earth_radius*cosZenith +  np.sqrt( (earth_radius*cosZenith)**2 - earth_radius**2 + (earth_radius+production_height+detector_depth)**2 )

        Osz = threeFlavorVacuumOsc()
        iOsc = 0
        if np.abs(pdg) == 14:
            iOsc = 1
        if np.abs(pdg) == 16:
            iOsc = 2

        flux_weight = 0.
        if pdg>0.:
            for f,i in zip(['nuE','nuMu'],[0,1]):
                flux_weight += Honda.getFlux(f,iCosZenithBin-1)[1][iEnergyBin-1] * Osz.posc(i, iOsc, baseline, E)

        if pdg<0.:
            for f,i in zip(['nuE','nuMu'],[0,1]):
                flux_weight += Honda.getFlux(f+'bar',iCosZenithBin-1)[1][iEnergyBin-1] * Osz.posc(i, iOsc, baseline, E)

        return flux_weight

    getHondaWeight = np.vectorize(getHondaWeight)

    print '...done'

    #######################
    # look at data
    #######################

    print 'process data'

    simWeight = np.array(input_data['weight'], dtype=float)

    atmosflux_weight = np.array(getHondaWeight(input_data['true_energy'],np.cos(input_data['true_zenith']),input_data['pdg']),dtype=float)
    weight = simWeight * atmosflux_weight

    oversampled_weight, oversampled_pdg, oversampled_energy_true, oversampled_energy_reco, oversampled_psi_true, oversampled_psi_reco, oversampled_RA_reco, oversampled_Dec_reco = oversample(weight, np.array(input_data['pdg']),input_data['true_energy'], input_data['reco_energy'], input_data['true_zenith'], input_data['reco_zenith'], input_data['true_azimuth'], input_data['reco_azimuth'], nOversampling)

    totalWeight = oversampled_weight
    print '...done'
     
    
print 'create PDF with binning', binning
    
binDef = {}
binDef['RA-DEC'] = [np.linspace(0,2*np.pi,360+1),np.linspace(-np.pi/2.,np.pi/2.,180+1)]
binDef['Psi-E'] = [np.linspace(0,np.pi,180+1),np.logspace(-1,3,4*20+1)]    
    
# make sure RA is in the right scale:
oversampled_RA_reco[oversampled_RA_reco<0.] = oversampled_RA_reco[oversampled_RA_reco<0.] + 2.*np.pi
    
if binning=='RA-DEC':
    pdf,xedges,yedges = np.histogram2d(oversampled_RA_reco, oversampled_Dec_reco,
               bins = binDef[binning],
               weights= totalWeight
    )
    pdf_quad,xedges,yedges = np.histogram2d(oversampled_RA_reco, oversampled_Dec_reco,
               bins = binDef[binning],
               weights= totalWeight**2
    )

if binning=='Psi-E':
    pdf,xedges,yedges = np.histogram2d(oversampled_psi_reco, oversampled_energy_reco,
               bins = binDef[binning],
               weights= totalWeight
    )
    pdf_quad,xedges,yedges = np.histogram2d(oversampled_psi_reco, oversampled_energy_reco,
               bins = binDef[binning],
               weights= totalWeight**2
    )

print '...done'

export = {}
export['pdf']=pdf
export['pdf_quad']=pdf_quad
export['bins_x'] = xedges
export['bins_y'] = yedges

PickleOutput = open(out_file, 'wb')
pickle.dump(export, PickleOutput)
PickleOutput.close()

print 'All done and PDF saved!'
