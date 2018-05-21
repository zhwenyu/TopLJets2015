#!/bin/bash

vars=(chmult chflux chfluxz chavgpt chavgpz chrecoil detST sphericity sphericity_2 aplanarity aplanarity_2 C C_2 D D_2)
vars=(chmult chflux chfluxz chavgpt chavgpz chrecoil sphericity aplanarity C D)
systs=(`cat store/TOP-17-015/C/inc/unfold/mean_summary.dat | awk '{print $1}'`)
header="Systs"
for i in ${vars[@]}; do 
    header="${header} & ${i}"
done
echo $header
for s in ${systs[@]}; do
    svals=${s}

    ncol=3
    if [ "$s" == "Trk." ]; then
        ncol=4
    fi
        
    for i in ${vars[@]}; do 
        val=`grep ${s} store/TOP-17-015/${i}/inc/unfold/mean_summary.dat | awk -v ncol=${ncol} '{print $ncol}'`; 
        svals="${svals} & ${val}";
    done
    echo $svals
done
