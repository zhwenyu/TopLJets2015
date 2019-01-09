#!/bin/bash

FITSCRIPT=test/analysis/top17010/createFit.py
OUTDIR=fit_results

function run_fit {
    echo "Running fit for $1"
    cd $1
    sh runFit.sh
    cd -
}


$FITSCRIPT -o $OUTDIR/em_inc -a \
    em_datacards/nom/datacard.dat 
run_fit $OUTDIR/em_inc

$FITSCRIPT -o $OUTDIR/em_ptcats -a \
    emhighpt=emhighpt_datacards/nom/datacard.dat emlowpt=emlowpt_datacards/nom/datacard.dat 
run_fit $OUTDIR/em_ptcats

$FITSCRIPT -o $OUTDIR/em_ptbcats -a \
    emhighpt2b=emhighpt2b_datacards/nom/datacard.dat emhighpt1b=emhighpt1b_datacards/nom/datacard.dat \
    emlowpt2b=emlowpt2b_datacards/nom/datacard.dat emlowpt1b=emlowpt1b_datacards/nom/datacard.dat 
run_fit $OUTDIR/em_ptbcats

$FITSCRIPT -o $OUTDIR/2l_inc -a \
    em=em_datacards/nom/datacard.dat ee=ee_datacards/nom/datacard.dat mm=mm_datacards/nom/datacard.dat 
run_fit $OUTDIR/2l_inc