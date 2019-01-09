#!/bin/bash

FITSCRIPT=test/analysis/top17010/createFit.py
OUTDIR=store/TOP17010/fit_results

function run_fit {
    echo "Running fit for $1"
    cd $1
    sh runFit.sh
    cd -
}

#loop over all the available scenarios
a=(`ls em_datacards`)
for i in ${a[@]}; do 

    baseOpts=""
    if [ ${i} == "nom" ]; then
        baseOpts="${baseOpts} -t 100 -a";
    fi
 
    #emu inclusive
    out=$OUTDIR/em_inc/${i};
    $FITSCRIPT -o ${out} ${baseOpts} em_datacards/${i}/datacard.dat;
    run_fit ${out};

    #emu+ee+mm inclusive
    out=$OUTDIR/dil_inc/${i};
    $FITSCRIPT -o ${out} ${baseOpts} \
        em=em_datacards/${i}/datacard.dat ee=ee_datacards/${i}/datacard.dat mm=mm_datacards/${i}/datacard.dat;
    run_fit ${out};

done





#$FITSCRIPT -o $OUTDIR/em_ptcats -a \
#    emhighpt=emhighpt_datacards/nom/datacard.dat emlowpt=emlowpt_datacards/nom/datacard.dat 
#run_fit $OUTDIR/em_ptcats

#$FITSCRIPT -o $OUTDIR/em_ptbcats -a \
#    emhighpt2b=emhighpt2b_datacards/nom/datacard.dat emhighpt1b=emhighpt1b_datacards/nom/datacard.dat \
#    emlowpt2b=emlowpt2b_datacards/nom/datacard.dat emlowpt1b=emlowpt1b_datacards/nom/datacard.dat 
#run_fit $OUTDIR/em_ptbcats

