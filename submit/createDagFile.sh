#!/bin/sh

masses=(1 2 4 6 8 10 20 40 60 80 100 150 200)

channels=('tau' 'b' 'nuMu' 'mu')
profiles=('NFW' 'Burkert')
binning=('Psi-E')

for b in "${binning[@]}"; do

    JOBIDBGPDF=background-PDF-$b
    echo JOB $JOBIDBGPDF PDF.submit
    echo VARS $JOBIDBGPDF JOBNAME=\"$JOBIDBGPDF\" TYPE=\"background\" CHANNEL=\"none\" PROFILE=\"none\" MASS=\"none\" OVERSAMPLING=\"100\" BINNING=\"$b\"

    for p in "${profiles[@]}"; do
        for c in "${channels[@]}"; do
            for m in "${masses[@]}"; do
                JOBIDPDF=$c-$m-GeV-PDF-$b-$p
                echo JOB $JOBIDPDF PDF.submit
                echo VARS $JOBIDPDF JOBNAME=\"$JOBIDPDF\" TYPE=\"signal\" CHANNEL=\"$c\" PROFILE=\"$p\" MASS=\"$m\" OVERSAMPLING=\"100\" BINNING=\"$b\"

                JOBIDLLH=$c-$m-GeV-LLH-$b-$p
                echo JOB $JOBIDLLH LLH.submit
                echo VARS $JOBIDLLH JOBNAME=\"$JOBIDLLH\" CHANNEL=\"$c\" PROFILE=\"$p\" MASS=\"$m\" OVERSAMPLING=\"100\" BINNING=\"$b\" REBINX=\"5\" REBINY=\"2\" CONF=\"90\"
               
                echo PARENT $JOBIDBGPDF $JOBIDPDF CHILD $JOBIDLLH

            done
        done
    done
done


c='tau'
b='RA-DEC'

JOBIDBGPDF=background-PDF-$b
echo JOB $JOBIDBGPDF PDF.submit
echo VARS $JOBIDBGPDF JOBNAME=\"$JOBIDBGPDF\" TYPE=\"background\" CHANNEL=\"none\" PROFILE=\"none\" MASS=\"none\" OVERSAMPLING=\"100\" BINNING=\"$b\"

for p in "${profiles[@]}"; do
    for m in "${masses[@]}"; do
        JOBIDPDF=$c-$m-GeV-PDF-$b-$p
        echo JOB $JOBIDPDF PDF.submit
        echo VARS $JOBIDPDF JOBNAME=\"$JOBIDPDF\" TYPE=\"signal\" CHANNEL=\"$c\" PROFILE=\"$p\" MASS=\"$m\" OVERSAMPLING=\"100\" BINNING=\"$b\"
        
        JOBIDLLH=$c-$m-GeV-LLH-$b-$p
        echo JOB $JOBIDLLH LLH.submit
        echo VARS $JOBIDLLH JOBNAME=\"$JOBIDLLH\" CHANNEL=\"$c\" PROFILE=\"$p\" MASS=\"$m\" OVERSAMPLING=\"100\" BINNING=\"$b\" REBINX=\"36\" REBINY=\"18\" CONF=\"90\"
        
        echo PARENT $JOBIDBGPDF $JOBIDPDF CHILD $JOBIDLLH
        
        
    done
done
