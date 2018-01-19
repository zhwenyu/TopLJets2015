#!/bin/bash

echo "Usage: ./getEMBuildupPlots.sh <STEP> <DIR>"
echo ""

STEP=$1
DIR=$2
queue=2nd

who=`whoami`
outdir=/afs/cern.ch/work/${who:0:1}/${who}/TOP-17-010-final/${DIR}/

wid=(20 40 50 60 70 80 90 
     100 110 120 130 140 150 
     160 180 200 220 240 260 
     280 300 350 400)        
#wid=(160)
lbCat=(highpt lowpt)
lfs=(EE EM MM)
cat=(1b 2b)
dists=(incmlb)

# Helpers: getting comma-separated lists of variables
function join { local IFS=','; echo "$*"; }
distStr="$(join ${dists[@]})"
lbcStr="$(join ${lbCat[@]})"
lfsStr="$(join ${lfs[@]})"
catStr="$(join ${cat[@]})"
widStr="$(join ${wid[@]})"

# CMSSW locations
CMSSW_7_4_7dir=~/CMSSW_7_4_7/src/
CMSSW_7_6_3dir=~/CMSSW_8_0_26_patch1/src/

# highptEM2b 
# --> +++++ highptEM1b
# --> +++++  lowptEM2b +  lowptEM1b
# --> +++++ highptEE2b + highptEE1b + lowptEE2b + lowptEE1b  
# --> +++++ highptMM2b + highptMM1b + lowptMM2b + lowptMM1b  

for dist in ${dists[*]}; do

    mkdir -p ${outdir}/EMBuildup_${dist}

for twid in ${wid[*]}; do
    cd ${CMSSW_7_4_7dir}
    eval `scramv1 runtime -sh`

    cardsdir=${outdir}/hypotest_100vs${twid}_100pseudodata/
    extdir=${outdir}/EMBuildup_${dist}/altWid_w${twid}
    mkdir -p ${extdir}


    echo "Starting buildup for ${dist} ${twid}"

    case $STEP in 
        CARDS)
        echo " - merging cards..."
            # high2b EM
            python ${CMSSW_7_4_7dir}/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
                EM2bhighpt=${cardsdir}/datacard_EM2bhighpt.dat \
                > ${extdir}/datacard_step1.dat

            # high1b + high2b EM
            python ${CMSSW_7_4_7dir}/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
                EM2bhighpt=${cardsdir}/datacard_EM2bhighpt.dat \
                EM1bhighpt=${cardsdir}/datacard_EM1bhighpt.dat \
                > ${extdir}/datacard_step2.dat

            # full EM ch
            python ${CMSSW_7_4_7dir}/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
                EM2bhighpt=${cardsdir}/datacard_EM2bhighpt.dat \
                EM1bhighpt=${cardsdir}/datacard_EM1bhighpt.dat \
                 EM2blowpt=${cardsdir}/datacard_EM2blowpt.dat \
                 EM1blowpt=${cardsdir}/datacard_EM1blowpt.dat \
                > ${extdir}/datacard_step3.dat

            # 2 ch of 3
            python ${CMSSW_7_4_7dir}/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
                EM2bhighpt=${cardsdir}/datacard_EM2bhighpt.dat \
                EM1bhighpt=${cardsdir}/datacard_EM1bhighpt.dat \
                 EM2blowpt=${cardsdir}/datacard_EM2blowpt.dat \
                 EM1blowpt=${cardsdir}/datacard_EM1blowpt.dat \
                EE2bhighpt=${cardsdir}/datacard_EE2bhighpt.dat \
                EE1bhighpt=${cardsdir}/datacard_EE1bhighpt.dat \
                 EE2blowpt=${cardsdir}/datacard_EE2blowpt.dat \
                 EE1blowpt=${cardsdir}/datacard_EE1blowpt.dat \
                > ${extdir}/datacard_step4.dat

            # all ch
            python ${CMSSW_7_4_7dir}/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
                EM2bhighpt=${cardsdir}/datacard_EM2bhighpt.dat \
                EM1bhighpt=${cardsdir}/datacard_EM1bhighpt.dat \
                 EM2blowpt=${cardsdir}/datacard_EM2blowpt.dat \
                 EM1blowpt=${cardsdir}/datacard_EM1blowpt.dat \
                EE2bhighpt=${cardsdir}/datacard_EE2bhighpt.dat \
                EE1bhighpt=${cardsdir}/datacard_EE1bhighpt.dat \
                 EE2blowpt=${cardsdir}/datacard_EE2blowpt.dat \
                 EE1blowpt=${cardsdir}/datacard_EE1blowpt.dat \
                MM2bhighpt=${cardsdir}/datacard_MM2bhighpt.dat \
                MM1bhighpt=${cardsdir}/datacard_MM1bhighpt.dat \
                 MM2blowpt=${cardsdir}/datacard_MM2blowpt.dat \
                 MM1blowpt=${cardsdir}/datacard_MM1blowpt.dat \
                > ${extdir}/datacard_step5.dat


            echo " - making workspaces..."

            for i in 1 2 3 4 5
            do
                cmd="text2workspace.py ${extdir}/datacard_step${i}.dat -P"
                cmd="${cmd} HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest"
                cmd="${cmd} -m 172.5 --PO verbose --PO altSignal=w${twid} --PO muFloating"
                cmd="${cmd} -o ${extdir}/workspace_step${i}.root"

                bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${extdir}" "${cmd}" 
            done


        ;;
        SCANS)
            echo " - running scans"

            cd ${extdir}/
            for i in 1 2 3 4 5
            do
                for j in 0 1
                do
                    cmd="combine ${extdir}/workspace_step${i}.root -M MultiDimFit"
                    cmd="${cmd} -m 172.5 -P x --floatOtherPOI=1 --algo=grid --points=200"
                    cmd="${cmd} -t -1 --expectSignal=1 --setPhysicsModelParameters x=${j},r=1"
                    cmd="${cmd} -n scan_x${j}_step${i}"
                    
                    bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${extdir}" "${cmd}" 
                done
            done
        ;;
        STATS)
            echo "- getting likelihood scans"
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getLikelihoodScans.py \
                -i ${extdir}/ \
                -o ${extdir}/ \
                --lfs step1,step2,step3,step4,step5 \
                --wid ${twid} \
                --dist ${dist}
        ;;
    esac

done
done
