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
outdir=/afs/cern.ch/work/${who:0:1}/${who}/TOP-17-010-final-v2/
extdir=${outdir}/${EXT}
CMSSW_7_4_7dir=~/CMSSW_8_1_0/src/
CMSSW_7_6_3dir=~/CMSSW_8_0_26_patch1/src/

unblind=true
nPseudo=1000

#mas=(1725)
mas=(1710 1712 1714 1716 1718 1720 1722 1724 1725 1726 1728 1730 1732 1734 1736 1738 1740)

#wid=(100)        
wid=(20 40 50 60 70 80 90 100 110 120 130 140 150 160 180 200 220 240 260 280 300 350 400)        

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
massSt="$(join ${mas[@]})"

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
               "tWttInterf" "tWMEScale" "all") 
        #nuisanceGroups=("jes0" "jes1" "jes2" "jes3" "jes4" "jes5" "jes6" "jes7" "jes8" "jes9" "jes10" 
        #            "jes11" "jes12" "jes13" "jes14" "jes15" "jes16" "jes17" "jes18" "jes19" 
        #            "jes20" "jes21" "jes22" "jes23" "jes24" "jes25" "jes26" "jes27" "jes28")
        #cats=(
        #    "EM1blowpt"
        #    "EM1bhighpt"
        #    "EM2blowpt"
        #    "EM2bhighpt"
        #    "MM1blowpt"
        #    "MM1bhighpt"
        #    "MM2blowpt"
        #    "MM2bhighpt"
        #)

        echo "------------------------"
        python test/TopWidthAnalysis/getContour.py \
            --mass 172.5 -n Contour1D_${3} \
            -e _m1725vs1725_100pseudodata \
            -i ${extdir}/

        #for i in ${cats[*]} ; do
        #    echo ${i}
        #    python test/TopWidthAnalysis/getContour.py \
        #        --mass 172.5 -n Contour1D_${i} \
        #        -e _${3}pseudodata \
        #        -i ${extdir}__${i}/
        #    echo "------------------------"
        #done

        #for i in 40 80 120 160 200 240 280 ; do
        #    python test/TopWidthAnalysis/getContour.py \
        #        --mass 172.5 -n Contour1D_${i} \
        #        -e _${i}pseudodata \
        #        -i ${extdir}/
        #done

        #for i in 1 2 3 4 ; do
        #    python test/TopWidthAnalysis/getContour.py \
        #        --mass 172.5 -n Contour1D_main${3}_step${i} \
        #        -e _${3}pseudodata \
        #        -i ${extdir}_step${i}/
        #done

        #for nuisGroup in ${nuisanceGroups[@]} ; do
        #    echo "------------------------"
        #    echo "Limits for ${nuisGroup}:"
        #    python test/TopWidthAnalysis/getContour.py \
        #        --mass 172.5 -n Contour1D_${3}_${nuisGroup} \
        #        -e _${3}pseudodata \
        #        -i ${extdir}Frz_${nuisGroup}/
        #done
    ;;
########################### 2D LIMITS #####################################
    2D_LIMITS)
        python test/TopWidthAnalysis/getContour.py \
            --mass ${massSt} --wids ${widStr} -n Contour2D_$3 \
            -e "_100pseudodata" \
            -i ${extdir}/

        #nuisanceGroups=("lumi_13TeV" "DYnorm_*CH*"
        #       "ees" "mes" "jer" "ltag" "btag" "bfrag" "semilep"
        #       "pu" "tttoppt" "ttMEqcdscale" "ttPDF"
        #       "jes" "UE" "CR" 
        #       "hdamp" "ISR_tt" "FSR_tt" 
        #       "tWttInterf" "tWMEScale" ) #"all") 

        #externArr=(
        #    "RateMod lumi_13TeV,Up"
        #    "RateMod lumi_13TeV,Dn"
        #    "RateMod DYnorm_*CH*,Up,EE"
        #    "RateMod DYnorm_*CH*,Up,EM"
        #    "RateMod DYnorm_*CH*,Up,MM"
        #    "RateMod DYnorm_*CH*,Dn,EE"
        #    "RateMod DYnorm_*CH*,Dn,EM"
        #    "RateMod DYnorm_*CH*,Dn,MM"
        #    "RateMod Wnorm_th,Up"
        #    "RateMod Wnorm_th,Dn"
        #    "RateMod tWnorm_th,Up"
        #    "RateMod tWnorm_th,Dn"
        #    "RateMod VVnorm_th,Up"
        #    "RateMod VVnorm_th,Dn"
        #    "RateMod tbartVnorm_th,Up"
        #    "RateMod tbartVnorm_th,Dn"
        #    "FromSim UEdn,tbart" 
        #    "FromSim UEup,tbart" 
        #    "FromSim QCDbased,tbart" 
        #    "FromSim gluonmove,tbart" 
        #    "FromSim ERDon,tbart" 
        #    "FromSim hdampdn,tbart" 
        #    "FromSim hdampup,tbart" 
        #    "FromSim isrdn,tbart" 
        #    "FromSim isrup,tbart" 
        #    "FromSim isrdn,Singletop" 
        #    "FromSim isrup,Singletop" 
        #    "FromSim fsrdn,tbart" 
        #    "FromSim fsrup,tbart" 
        #    "FromSim fsrdn,Singletop" 
        #    "FromSim fsrup,Singletop" 
        #    "FromSim DS,Singletop" 
        #    "FromSim medn,Singletop" 
        #    "FromSim meup,Singletop" 
        #    "FromSim UEdn,tbart" 
        #    "FromSim UEup,tbart" 
        #)

	    #    for exterName in "${externArr[@]}"; do
        #        echo "------------------------"
        #        iNuisExt=${exterName//\*/}
        #        iNuisExt=${iNuisExt//,/_}
        #        iNuisExt=${iNuisExt// /}
        #        echo "Limits for ${iNuisExt}:"
        #        python test/TopWidthAnalysis/getContour.py \
        #            --mass ${massSt} --wids ${widStr} \
        #            -n Contour2D_${iNuisExt} \
        #            -e "" \
        #            -i ${extdir}_Frz_all_Ext_${iNuisExt}/
        #    done
    ;;
############################### PLOT_NUIS #################################
    PLOT_NUIS )
        cd ${CMSSW_7_4_7dir}
        eval `scramv1 runtime -sh` 
        cd - 

        nuisanceGroups=("sel,trig_*CH*" "lumi_13TeV" "DYnorm_*CH*" "Wnorm_th" 
               "tWnorm_th" "VVnorm_th" "tbartVnorm_th" 
               "ees" "mes" "jer" "ltag" "btag" "bfrag" "semilep"
               "pu" "tttoppt" "ttMEqcdscale" "ttPDF"
               "jes" "st_wid" "UE" "CR" 
               "hdamp" "ISR" "FSR" "mtop" 
               "tWttInterf" "tWMEScale" "all") 

        for ngrp in ${nuisanceGroups[*]} ; do
            iwidN=0
        for twid in ${wid[*]} ; do
            imassN=0
        for tmass in ${mas[*]} ; do
            nwid=${widNs[${iwidN}]}
            nmas=${masNs[${imassN}]}

            python test/TopWidthAnalysis/getNuisances.py \
                -i ${extdir}_only${ngrp}/hypotest_100vs${twid}_m1725vs${tmass}_${3}pseudodata/higgsCombineNuisancesRun.MultiDimFit.mH172.5.123456.root \
                -o ${extdir}_only${ngrp}/ \
                -n Nuisances_${twid}_${tmass}_incmlb_only${ngrp} \
                --extraText "#Gamma_{Alt.} = ${nwid} #times #Gamma_{SM} ; m_{Alt.} = ${nmas} GeV"

            let "imassN += 1"
        done
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
