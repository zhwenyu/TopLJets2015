#!/bin/bash

WHAT=$1; 
EXT=$2;
if [ "$#" -lt 3 ]; then 
    echo "steerTOPWidthAnalysis.sh <COMMAND> <DIRECTORY> <PSEUDODATA>";
    echo "  Possible <COMMAND>:";
    echo "        CLs             - calculates CLs using higgs_combine";
    echo "        TOYS            - plots toys and saves locally";
    echo "        QUANTILES       - plots the quantiles chart for all distributions, by lepton final state";
    echo "        PLOT_NUIS       - plots post-fit nuisances for all regimes"; 
    echo "        PLOT_SYSTS      - plots systematic variations for given <PSEUDODATA>"; 
    echo " ";
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=Y

queue=2nw
who=`whoami`
outdir=/afs/cern.ch/work/${who:0:1}/${who}/TOP-17-010-final/
extdir=${outdir}/${EXT}
CMSSW_7_4_7dir=~/CMSSW_7_4_7/src/
CMSSW_7_6_3dir=~/CMSSW_8_0_26_patch1/src/

unblind=true
nPseudo=1000

wid=(20 40 50 60 70 80 90 100 110 120 130 140 150 160 180 200 220 240 260 280 300 350 400)        
widNs=("0.20" "0.40" "0.50" "0.60" "0.70" "0.80" "0.90" 
       "1.00" "1.10" "1.20" "1.30" "1.40" "1.50" "1.60" "1.80"
       "2.00" "2.20" "2.40" "2.60" "2.80"
       "3.00" "3.50" "4.00") 

lbCat=(highpt lowpt)
lfs=(EE EM MM)
cat=(1b 2b)
dists=(incmlb)

RED='\e[31m'
NC='\e[0m'


# Helpers: getting comma-separated lists of variables
function join { local IFS=','; echo "$*"; }
distStr="$(join ${dists[@]})"
lbcStr="$(join ${lbCat[@]})"
lfsStr="$(join ${lfs[@]})"
catStr="$(join ${cat[@]})"
widStr="$(join ${wid[@]})"
nwidSt="$(join ${widNs[@]})"

mkdir -p $extdir

case $WHAT in
############################# LIMITS #####################################
    LIMITS )
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd -

        nuisanceGroups=("sel,trig_*CH*" "lumi_13TeV" "DYnorm_*CH*" "Wnorm_th" 
               "tWnorm_th" "VVnorm_th" "tbartVnorm_th" 
               "ees" "mes" "jer" "ltag" "btag" "bfrag" "semilep"
               "pu" "tttoppt" "ttMEqcdscale" "ttPDF"
               "jes" "st_wid" "UE" "CR" 
               "hdamp" "ISR" "FSR" "mtop" 
               "tWttInterf" "tWMEScale") 

        echo "------------------------"
        python test/TopWidthAnalysis/getContour.py \
            --mass 172.5 -n Contour1D_main \
            -i ${extdir}/

        #for i in 1 2 3 4 ; do
        #    python test/TopWidthAnalysis/getContour.py \
        #        --mass 172.5 -n Contour1D_main_step${i} \
        #        -i ${extdir}_step${i}/
        #done

        for nuisGroup in ${nuisanceGroups[@]} ; do
            echo "------------------------"
            echo "Limits for ${nuisGroup}:"
            python test/TopWidthAnalysis/getContour.py \
                --mass 172.5 -n Contour1D_${nuisGroup} \
                -i ${extdir}Frz_${nuisGroup}/
        done
    ;;
############################### CLs #######################################
    CLs ) # get CLs statistics from combine
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${extdir}
        for dist in ${dists[*]} ; do
        for twid in ${wid[*]} ; do

            # pre-fit 
            echo "Making CLs for ${twid} ${dist}"
            cmd="combine ${extdir}/hypotest_100vs${twid}_${3}pseudodata/workspace.root -M HybridNew --seed 8192 --saveHybridResult" 
            cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 8"
            cmd="${cmd} --clsAcc 0 --fullBToys --saveToys --saveWorkspace --generateExt=1 --generateNuis=0"
            cmd="${cmd} --expectedFromGrid 0.5 -n cls_prefit_exp"

            bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${extdir}/hypotest_100vs${twid}_${3}pseudodata/" "${cmd}" 


            # post-fit expected 
            echo "Making CLs for ${twid} ${dist}"
            cmd="combine ${extdir}/hypotest_100vs${twid}_${3}pseudodata/workspace.root -M HybridNew --seed 8192 --saveHybridResult" 
            cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 8"
            cmd="${cmd} --clsAcc 0 --fullBToys  --saveWorkspace --saveToys --frequentist"
            cmd="${cmd} --expectedFromGrid 0.5 -n cls_postfit_exp"

            bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${extdir}/hypotest_100vs${twid}_${3}pseudodata/" "${cmd}" 

            # post-fit observed
            if [[ ${unblind} == true ]] ; then 
                echo "Making CLs for ${twid} ${dist}"
                cmd="combine ${extdir}/hypotest_100vs${twid}_${3}pseudodata/workspace.root -M HybridNew --seed 8192 --saveHybridResult" 
                cmd="${cmd} -m 172.5  --testStat=TEV --singlePoint 1 -T ${nPseudo} -i 2 --fork 8"
                cmd="${cmd} --clsAcc 0 --fullBToys --saveWorkspace --saveToys --frequentist"
                cmd="${cmd} -n cls_postfit_obs"

                bsub -q ${queue} ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh "${extdir}/hypotest_100vs${twid}_${3}pseudodata/" "${cmd}" 
            fi
        done
        done
    ;;
############################### TOYS ######################################
    TOYS ) # get toys distributions from the pseudoexperiments
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh`
        cd ${extdir}
        for dist in ${dists[*]} ; do

            tWidI=0
        for twid in ${wid[*]} ; do
            widN=${widNs[${tWidI}]}

            # pre-fit expected 
            rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
            rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\("
            rootcmds="${rootcmds}172.5,1,\\\"x\\\",\\\"\\\",\\\"${widN}\\\",\\\"${dist}\\\",false,\\\"pre\\\"\)"

            cmd=""
            cmd="${cmd}root -l -q -b"
            cmd="${cmd} ${extdir}"
            cmd="${cmd}/hypotest_100vs${twid}_${3}pseudodata/higgsCombinecls_prefit_exp.HybridNew.mH172.5.8192.quant0.500.root"
            cmd="${cmd} ${rootcmds}"

            sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${extdir}/" "${cmd}"


            # post-fit expected 
            rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
            rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\("
            rootcmds="${rootcmds}172.5,1,\\\"x\\\",\\\"\\\",\\\"${widN}\\\",\\\"${dist}\\\",false,\\\"post\\\"\)"

            cmd=""
            cmd="${cmd}root -l -q -b"
            cmd="${cmd} ${extdir}"
            cmd="${cmd}/hypotest_100vs${twid}_${3}pseudodata/higgsCombinecls_postfit_exp.HybridNew.mH172.5.8192.quant0.500.root"
            cmd="${cmd} ${rootcmds}"

            sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                "${extdir}/hypotest_100vs${twid}_${3}pseudodata/" "${cmd}"

            # post-fit observed 
            if [[ ${unblind} == true ]] ; then 
                rootcmds="${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis"
                rootcmds="${rootcmds}/hypoTestResultTreeTopWid.cxx\("
                rootcmds="${rootcmds}172.5,1,\\\"x\\\",\\\"\\\",\\\"${widN}\\\",\\\"${dist}\\\",${unblind},\\\"obs\\\"\)"

                cmd=""
                cmd="${cmd}root -l -q -b"
                cmd="${cmd} ${extdir}"
                cmd="${cmd}/hypotest_100vs${twid}_${3}pseudodata/higgsCombinecls_postfit_obs.HybridNew.mH172.5.8192.root"
                cmd="${cmd} ${rootcmds}"

                sh ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/wrapPseudoexperiments.sh \
                        "${extdir}/hypotest_100vs${twid}_${3}pseudodata/" "${cmd}"
            fi

            cp ${extdir}/hypotest_100vs${twid}_${3}pseudodata/*obs*${twid}*${dist}*.{pdf,png} ${extdir}
            cp ${extdir}/hypotest_100vs${twid}_${3}pseudodata/*post*${twid}*${dist}*.{pdf,png} ${extdir}
            cp ${extdir}/hypotest_100vs${twid}_${3}pseudodata/*pre*${twid}*${dist}*.{pdf,png} ${extdir}
            cp ${extdir}/hypotest_100vs${twid}_${3}pseudodata/*obs*stats*.txt ${extdir}
            cp ${extdir}/hypotest_100vs${twid}_${3}pseudodata/*post*stats*.txt ${extdir}
            cp ${extdir}/hypotest_100vs${twid}_${3}pseudodata/*pre*stats*.txt ${extdir}

            let "tWidI += 1"
        done
        done
    ;;
    QUANTILES ) # plot quantiles distributions of all toys, get CLsPlot
            
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh` 

        cd ${extdir}
        rm statsPlots.root
        for dist in ${dists[*]} ; do

            ## Quantiles plot with pre-fit information 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getQuantilesPlot.py \
                -i ${extdir}/ -o ${extdir}/ \
                --wid ${widStr} \
                --nwid ${nwidSt} \
                --dist ${dist}  \
                --prep pre

            ## Quantiles plot with post-fit information 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getQuantilesPlot.py \
                -i ${extdir}/ -o ${extdir}/ \
                --wid ${widStr} \
                --nwid ${nwidSt} \
                --dist ${dist}  \
                --prep post \
                --unblind

            ## Quantiles plot with observed information 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getQuantilesPlot.py \
                -i ${extdir}/ -o ${extdir}/ \
                --wid ${widStr} \
                --nwid ${nwidSt} \
                --dist ${dist}  \
                --prep obs \
                --unblind
            
            # Get CLs plots for pre-fit expectations
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
                -i ${extdir}/ -o ${extdir}/ \
                --wid ${widStr} \
                --prep pre \
                --dist ${dist}
            
            # Get CLs plots for post-fit expectations 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
                -i ${extdir}/ -o ${extdir}/ \
                --wid ${widStr} \
                --dist ${dist} \
                --prep post \
                --addPre \
                --unblind
            
            # Get CLs plots for observed expectations 
            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getSeparationTables.py \
                -i ${extdir}/ -o ${extdir}/ \
                --wid ${widStr} \
                --dist ${dist} \
                --prep obs \
                --unblind

            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getCLsFromFit.py \
                -i ${extdir}/ \
                --dist ${dist} \
                --prep pre

            python ${CMSSW_7_6_3dir}/TopLJets2015/TopAnalysis/test/TopWidthAnalysis/getCLsFromFit.py \
                -i ${extdir}/ \
                --dist ${dist} \
                --prep post \
                --unblind
        done
    ;;
############################### PLOT_NUIS #################################
    PLOT_NUIS )
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh` 
        cd - 

        for dist in ${dists[*]} ; do
            iwidN=0
        for twid in ${wid[*]} ; do
            nwid=${widNs[${iwidN}]}

            python test/TopWidthAnalysis/getNuisances.py \
                -i ${extdir}/hypotest_100vs${twid}_${3}pseudodata/higgsCombinecls_prefit_exp.HybridNew.mH172.5.8192.quant0.500.root \
                -o ${extdir}/ \
                -n preNuisances_${twid}_${dist} \
                --extraText "#Gamma_{Alt.} = ${nwid} #times #Gamma_{SM}"

            python test/TopWidthAnalysis/getNuisances.py \
                -i ${extdir}/hypotest_100vs${twid}_${3}pseudodata/higgsCombinecls_postfit_exp.HybridNew.mH172.5.8192.quant0.500.root \
                -o ${extdir}/ \
                -n postNuisances_${twid}_${dist} \
                --extraText "#Gamma_{Alt.} = ${nwid} #times #Gamma_{SM}"

            python test/TopWidthAnalysis/getNuisances.py \
                -i ${extdir}/hypotest_100vs${twid}_${3}pseudodata/higgsCombinecls_postfit_obs.HybridNew.mH172.5.8192.root \
                -o ${extdir}/ \
                -n obsNuisances_${twid}_${dist} \
                --extraText "#Gamma_{Alt.} = ${nwid} #times #Gamma_{SM}"

            let "iwidN += 1"
        done
        done
    ;;
############################### PLOT_SYSTS #################################
    PLOT_SYSTS )
        catList="EE1bhighpt,EE2bhighpt,EE1blowpt,EE2blowpt,"
        catList="${catList}EM1bhighpt,EM2bhighpt,EM1blowpt,EM2blowpt,"
        catList="${catList}MM1bhighpt,MM2bhighpt,MM1blowpt,MM2blowpt"

        python test/TopWidthAnalysis/getShapeUncPlots.py \
            --input root://eoscms//eos/cms/store/cmst3/group/top/TOP-17-010-final/plotter/plotter.root \
            --tag exp \
            --wid $3 \
            --cats ${catList} \
            --proc "t#bar{t},Single top"
        python test/TopWidthAnalysis/getShapeUncPlots.py \
            --input root://eoscms//eos/cms/store/cmst3/group/top/TOP-17-010-final/plotter/plotter.root \
            --wid $3 \
            --tag gen \
            --cats ${catList} \
            --proc "t#bar{t}"
        python test/TopWidthAnalysis/getShapeUncPlots.py \
            --input root://eoscms//eos/cms/store/cmst3/group/top/TOP-17-010-final/plotter/plotter.root \
            --systInput root://eoscms//eos/cms/store/cmst3/group/top/TOP-17-010-final/plotter/syst_plotter.root \
            --wid $3 \
            --tag sys \
            --cats ${catList} \
            --proc "t#bar{t},Single top"

    ;;
esac
