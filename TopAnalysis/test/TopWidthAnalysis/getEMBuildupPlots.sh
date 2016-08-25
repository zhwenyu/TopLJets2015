#!/bin/bash

ERA=$1
STEP=$2
queue=2nd

outdir=/afs/cern.ch/work/e/ecoleman/public/TopWidth/TopWidth_${ERA}_808p1
cardsdir=${outdir}/datacards

lumi=11377
case $ERA in
    era2015)
    lumi=2267.84
    ;;
esac

lfs=(EE EM MM)
wid=(0p5w 1p0w 1p5w 2p0w 2p5w 3p0w 3p5w 4p0w 4p5w 5p0w)
dists=(incmlb sncmlb mt2mlb)
cat=(1b 2b)
lbCat=(highpt lowpt)
unblind=false
nPseudo=1000

RED='\e[31m'
NC='\e[0m'

CMSSW_7_4_7dir=~/CMSSW_7_4_7/src/
CMSSW_7_6_3dir=~/CMSSW_8_0_8_patch1_2/src/

# highptEM2b 
# --> +++++ highptEM1b
# --> +++++ lowptEM2b + lowptEM1b
# --> +++++ highptEE2b + highptEE1b + lowptEE2b + lowptEE1b  
# --> +++++ highptMM2b + highptMM1b + lowptMM2b + lowptMM1b  

for dist in ${dists[*]}; do

    mkdir ${outdir}/EMsteps_${dist}

for twid in ${wid[*]}; do
    cd ${CMSSW_7_4_7dir}
    eval `scramv1 runtime -sh`

    echo "Starting buildup for ${dist} ${twid}"

    # TODO: Add in 1p0w hypo tests 
    if [[ "${twid}" == "1p0w" ]] ; then
        continue
    fi

    case $STEP in 
        CARDS)
        echo " - merging cards..."
            # high2b EM
            python ${CMSSW_7_4_7dir}/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
                highptEM2b=${cardsdir}/datacard__${twid}_highptEM2b_${dist}.dat \
                > ${cardsdir}/datacard__${twid}_${dist}_step1.dat

            # high1b + high2b EM
            python ${CMSSW_7_4_7dir}/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
                highptEM2b=${cardsdir}/datacard__${twid}_highptEM2b_${dist}.dat \
                highptEM1b=${cardsdir}/datacard__${twid}_highptEM1b_${dist}.dat \
                > ${cardsdir}/datacard__${twid}_${dist}_step2.dat


            # full EM ch
            python ${CMSSW_7_4_7dir}/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
                highptEM2b=${cardsdir}/datacard__${twid}_highptEM2b_${dist}.dat \
                highptEM1b=${cardsdir}/datacard__${twid}_highptEM1b_${dist}.dat \
                lowptEM2b=${cardsdir}/datacard__${twid}_lowptEM2b_${dist}.dat \
                lowptEM1b=${cardsdir}/datacard__${twid}_lowptEM1b_${dist}.dat \
                > ${cardsdir}/datacard__${twid}_${dist}_step3.dat

            # 2 ch of 3
            python ${CMSSW_7_4_7dir}/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
                highptEM2b=${cardsdir}/datacard__${twid}_highptEM2b_${dist}.dat \
                highptEM1b=${cardsdir}/datacard__${twid}_highptEM1b_${dist}.dat \
                 lowptEM2b=${cardsdir}/datacard__${twid}_lowptEM2b_${dist}.dat \
                 lowptEM1b=${cardsdir}/datacard__${twid}_lowptEM1b_${dist}.dat \
                highptEE2b=${cardsdir}/datacard__${twid}_highptEE2b_${dist}.dat \
                highptEE1b=${cardsdir}/datacard__${twid}_highptEE1b_${dist}.dat \
                 lowptEE2b=${cardsdir}/datacard__${twid}_lowptEE2b_${dist}.dat \
                 lowptEE1b=${cardsdir}/datacard__${twid}_lowptEE1b_${dist}.dat \
                > ${cardsdir}/datacard__${twid}_${dist}_step4.dat

            # all ch
            python ${CMSSW_7_4_7dir}/HiggsAnalysis/CombinedLimit/scripts/combineCards.py \
                highptEM2b=${cardsdir}/datacard__${twid}_highptEM2b_${dist}.dat \
                highptEM1b=${cardsdir}/datacard__${twid}_highptEM1b_${dist}.dat \
                 lowptEM2b=${cardsdir}/datacard__${twid}_lowptEM2b_${dist}.dat \
                 lowptEM1b=${cardsdir}/datacard__${twid}_lowptEM1b_${dist}.dat \
                highptEE2b=${cardsdir}/datacard__${twid}_highptEE2b_${dist}.dat \
                highptEE1b=${cardsdir}/datacard__${twid}_highptEE1b_${dist}.dat \
                 lowptEE2b=${cardsdir}/datacard__${twid}_lowptEE2b_${dist}.dat \
                 lowptEE1b=${cardsdir}/datacard__${twid}_lowptEE1b_${dist}.dat \
                highptMM2b=${cardsdir}/datacard__${twid}_highptMM2b_${dist}.dat \
                highptMM1b=${cardsdir}/datacard__${twid}_highptMM1b_${dist}.dat \
                 lowptMM2b=${cardsdir}/datacard__${twid}_lowptMM2b_${dist}.dat \
                 lowptMM1b=${cardsdir}/datacard__${twid}_lowptMM1b_${dist}.dat \
                > ${cardsdir}/datacard__${twid}_${dist}_step5.dat


            echo " - making workspaces..."

            for i in 1 2 3 4 5
            do
                text2workspace.py ${cardsdir}/datacard__${twid}_${dist}_step${i}.dat -P \
                    HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest \
                    -m 172.5 --PO verbose --PO altSignal=${twid} --PO muFloating \
                    -o ${outdir}/EMsteps_${dist}/${twid}_step${i}.root
            done


        ;;
        SCANS)
            echo " - running scans"

            cd ${outdir}/EMsteps_${dist}/
            for i in 1 2 3 4 5
            do
                for j in 0 1
                do
                    cmd="combine ${outdir}/EMsteps_${dist}/${twid}_step${i}.root -M MultiDimFit"
                    cmd="${cmd} -m 172.5 -P x --floatOtherPOI=1 --algo=grid --points=200"
                    cmd="${cmd} -t -1 --expectSignal=1 --setPhysicsModelParameters x=${j},r=1"
                    cmd="${cmd} -n x${j}_scan_Asimov_${twid}_step${i}"
                    
                    bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${outdir}/EMsteps_${dist}/" "${cmd}" 
                done
            done
        ;;
        CLs)
            echo " - producing CLs stats"

            cd ${outdir}/EMsteps_${dist}/
            for i in 1 2 3 4 5
            do
                cmd="combine ${outdir}/EMsteps_${dist}/${twid}_step${i}.root -M HybridNew --seed 8192 --saveHybridResult"
                cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 6"
                cmd="${cmd} --clsAcc 0 --fullBToys  --generateExt=1 --generateNuis=0" 
                cmd="${cmd} --expectedFromGrid 0.5 -n x_pre-fit_exp__${twid}_step${i}" 
                #cmd="${cmd} &> ${outdir}/EMsteps_${dist}/x_pre-fit_exp__${twid}_step${i}.log"
                
                bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${outdir}/EMsteps_${dist}/" "${cmd}" 
            done
        ;;
        TOYS)
            echo " - getting toys"

            cd ${outdir}/EMsteps_${dist}/
            for i in 1 2 3 4 5
            do
                rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
                rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx(\"x_pre-fit_exp__${twid}_step${i}.qvals.root"
                rootcmds="${rootcmds}\",172.5,1,\"x\",${nPseudo},\"step${i}\",\"${twid}\",\"${dist}\",${unblind})"
                # All datacards
                root -l -q -b higgsCombinex_pre-fit_exp__${twid}_step${i}.HybridNew.mH172.5.8192.quant0.500.root \
                    ${rootcmds}
            done

            echo " - getting quantiles plot"

            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getQuantilesPlot.py \
                -i ${outdir}/EMsteps_${dist}/ \
                -o ${outdir}/EMsteps_${dist}/${twid}_ \
                --lfs step1,step2,step3,step4,step5 \
                --dist ${dist} \
                --wid ${twid}   \
                --labelWidth \
                --axisOverwrite "(1),(2),(3),(4),(5)"
        ;;
    esac

done

    # get likelihood scan plots
    case $STEP in
        STATS)
            echo "- getting likelihood scans"
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getLikelihoodScans.py -i ${outdir}/EMsteps_${dist}/ -o ${outdir}/EMsteps_${dist}/
        ;;
    esac

done
