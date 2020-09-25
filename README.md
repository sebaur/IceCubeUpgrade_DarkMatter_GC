# IceCube Upgrade: sensitivity for low mass dark kmatter in the Galactic Halo
Analysis code for the search of annihilating dark matter in the center of the Galaxy with the upgraded IceCube detector

Analysis location: `/data/user/sbaur/projects/finalized/Upgrade_DarkMatter/`

## Code/File structure:
- README.md : this file
- env.sh : file to define environment variables needed by the scripts, needs to be modified according to the place this analysis is copied to
- PDFs : folder containing all PDFs in pickle format
- sensitivity : folder containing calculated sensitivities in pickle format
- resources : folder containing code and data files needed for the calculation, e.g. DM halo profiles, the PPPC4 and Honda2015 tables
- python : various python code for the analysis
  - LLHAnalyser.py: Analyser class for the likelihood analysis
  - createPDFs.py : code for the calculation of PDFs with parameters
    - `-t type` (background/signal)
    - `-o oversampling_factor`
    - `-c channel` annihilation channel (signal only)
    - `-m mass` (signal only)
    - `-p profile` DM profile NFW/Burkert (signal only)
    - `-b binning` how to bin the PDFS, i.e. in RA-DEC or Psi-E  
    
  - Sensitivity.py : code for sensitivity calculation with parameters:
    - `-c channel` annihilation/decay channel
    - `-m mass`
    - `-p profile` DM profile NFW/Burkert
    - `-t livetime` analysis livetime in years
    - `-b binning` chosen binning (RA-DEC or Psi-E)
    - `-o oversampling_factor`
    - `-n confidence_level`
    - `-l likelihood_method` (Poisson or Effective)
    - `-x x_rebinning`
    - `-y y_rebinning`

- plotting : jupyter notebooks with examples for plotting PDFs or sensitivities
- submit : example scripts for condor/dagman submission


## Prerequisites
- icerec build, here from /data/user/sbaur/metaprojects/icerec/build/
- python 2, from /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1
- python modules: iminuit, scipy, astropy, OptionParser

## Workflow / examples

1. PDF creation for signal and background, specific for each channel and mass; example for 50 GeV annihilation to tau pairs:
```
./createPDFs.py -t signal -p NFW -c tau -m 50 -o 100 -b Psi-E
./createPDFs.py -t background -o 100 -b Psi-E
```

2. Calculate sensitivity for each channel and mass; example for 50 GeV annihilation to tau pairs:
```
./Sensitivity.py -c tau -m 50 -p NFW -t 3 -b Psi-E -o 100 -n 90 -l Poisson -x 5 -y 2
```
