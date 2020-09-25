#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /data/user/sbaur/metaprojects/icerec/build/

import numpy as np
import scipy
import pickle
import sys, os
import shutil
from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-c", "--channel",default="nu",
                  dest="CHANNEL", help="Annihilation channel")
parser.add_option("-p", "--profile",default="NFW",
                  dest="PROFILE", help="DM Profile")
parser.add_option("-m", "--mass", default='0',
                  dest="MASS", help="DM mass")
parser.add_option("-t", "--livetime", default='3',
                  dest="TIME", help="livetime")
parser.add_option("-o", "--oversampling", default='100',
                  dest="OVERSAMPLING", help="Oversampling")
parser.add_option("-b", "--binning", default='RA-DEC',
                  dest="BINNING", help="Bins for PDF")
parser.add_option("-l", "--llh", default='Poisson',
                  dest="LLH", help="likelihood type")
parser.add_option("-x", "--rebinx", default='1',
                  dest="REBINX", help="rebinning x")
parser.add_option("-y", "--rebiny", default='1',
                  dest="REBINY", help="rebinning y")

parser.add_option("-n", "--confidence", default='90',
                  dest="CONFIDENCE", help="confidence level")


(options,args) = parser.parse_args()
        
channel = options.CHANNEL
profile = options.PROFILE
binning = options.BINNING
mass = float(options.MASS)
livetime = float(options.TIME)*3600.*24.*356.
nOversampling = int(options.OVERSAMPLING)
rebin_x = int(options.REBINX)
rebin_y = int(options.REBINY)
LLH_type = options.LLH


source = 'source /data/user/sbaur/projects/PrepareForGit/Upgrade_DarkMatter/env.sh'
dump = '/usr/bin/python -c "import os,pickle;print pickle.dumps(os.environ)"'
penv = os.popen('%s && %s' %(source,dump))
env = pickle.loads(penv.read())
os.environ = env

base_path = os.environ['ANALYSIS_BASE_PATH']
PDF_path = base_path+'/PDFs/'

out_filename = 'Sensitivity_'+LLH_type+'_'+channel+'_'+str(int(mass))+'_'+profile+'_oversampling'+str(nOversampling)+'_'+binning+'_rebin'+str(rebin_x)+'-'+str(rebin_y)+'.npy'

print 'save result to file', base_path+out_filename
if os.path.isfile(base_path+'/sensitivity/'+out_filename):
    print " File exists already -- won't do anything!"
    sys.exit(0)

sys.path.append(base_path+'/python/')
import LLHAnalyser


def merge_bins(hist, binsmerge_x, binsmerge_y):
    '''
    This function merges the bins in a histogram
    '''
    
    assert np.shape(hist)[0]%binsmerge_x ==0 , 'merging not possible'
    assert np.shape(hist)[1]%binsmerge_y ==0 , 'merging not possible'
    
    len_new_x = int(np.shape(hist)[0] / binsmerge_x)
    len_new_y = int(np.shape(hist)[1] / binsmerge_y)
    
    new_hist = np.zeros((len_new_x,len_new_y))
    
    for j in range(binsmerge_x):
        for k in range(len_new_x):
            for j2 in range(binsmerge_y):
                for k2 in range(len_new_y):
                    new_hist[k,k2] += hist[j+k*binsmerge_x , j2+k2*binsmerge_y]
                    
    return new_hist


bg_file = 'Background_PDF_Honda2015_oversampling'+str(int(nOversampling))+'_'+binning+'.pkl'
with open(PDF_path+'/background/'+bg_file, 'r') as bg_pkl_file:
    bgFileData = pickle.load(bg_pkl_file)  
    bgpdf = merge_bins(bgFileData['pdf'],rebin_x ,rebin_y)
    bgpdf_unc = merge_bins(bgFileData['pdf_quad'],rebin_x ,rebin_y)


sig_file = 'Signal_PDF_'+profile+'_'+str(int(mass))+'GeV_'+channel+'_oversampling'+str(int(nOversampling))+'_'+binning+'.pkl'
with open(PDF_path+'/signal/'+sig_file, 'r') as sig_pkl_file:
    sigFileData = pickle.load(sig_pkl_file)
    sigpdf = merge_bins(sigFileData['pdf'], rebin_x, rebin_y)
    sigpdf_unc = merge_bins(sigFileData['pdf_quad'], rebin_x, rebin_y)

       
# init LLH analysis    
analysis = LLHAnalyser.Profile_Analyser()

analysis.loadBackgroundPDF(bgpdf.flatten()*livetime)
analysis.loadSignalPDF(sigpdf.flatten()*livetime)

if LLH_type == 'Effective':
    analysis.loadUncertaintyPDFs(bgpdf_unc.flatten()*(livetime**2),sigpdf_unc.flatten()*(livetime**2))

print 'Using',LLH_type,'likelihood method.'
analysis.setLLHtype(LLH_type)

Ntrials = 10000
sens = analysis.CalculateSensitivity(Ntrials, int(options.CONFIDENCE))

print 'Median sensitivity signal fraction',sens['median']

xs_norm = np.sum(bgpdf)/np.sum(sigpdf)
print 'xs_norm=',xs_norm

sens['mass'] = mass
sens['error_68_low'] = sens['error_68_low']*xs_norm
sens['error_68_high'] = sens['error_68_high']*xs_norm
sens['error_95_low'] = sens['error_95_low']*xs_norm
sens['error_95_high'] = sens['error_95_high']*xs_norm
sens['median'] = sens['median']*xs_norm

print 'Saving output file'

np.save(base_path+'/sensitivity/'+out_filename,sens)

print '...done'

