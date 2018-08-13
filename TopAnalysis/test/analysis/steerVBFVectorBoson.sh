#!/bin/bash

WHAT=$1; 
EXTRA=$2
if [ "$#" -lt 1 ]; then 
    echo "steerVBFVectorBoson.sh <SEL/MERGE/PLOT/WWW> [extra]";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGE        - merge output (if given \"extra\" is appended to the directory)"
    echo "        PLOT         - make plots (if given \"extra\" is appended to the directory)"
    echo "        WWW          - move plots to web-based (if given \"extra\" is appended to the directory)"
    exit 1; 
fi

#to run locally use local as queue + can add "--njobs 8" to use 8 parallel jobs
queue=workday
githash=fbc74ae
eosdir=/store/cmst3/group/top/RunIIFall17/${githash}
githash2018=4ad3a45
eosdir2018=/store/cmst3/group/top/RunIISpring18/${githash2018}
githashJetHT=dc6c6d0
eosdirJetHT=/store/cmst3/user/ajafari/RunIIFall17/${githashJetHT}
fulllumi=41367
vbflumi=7661
fulllumi2018=5300
lumiUnc=0.025
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBFVectorBoson
wwwdir=~/www/VBFVectorBoson


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
        # input=${eosdir}/MC13TeV_DY50toInf/MergedMiniEvents_0_ext0.root
        # output=MC13TeV_DY4Jets50toInf.root
        # tag="--tag MC13TeV_DY50toInf"

        input=${eosdirJetHT}/Data13TeV_JetHT_2017D/MergedMiniEvents_4_ext0.root #/MC13TeV_AJJ_EWK_INT_LO_mjj500_dr04/MC13TeV_AJJ_EWK_INT_LO_mjj500_dr04.root
        output=Data13TeV_JetHT_2017D.root #MC13TeV_AJJ_EWK_INT_LO_mjj500_dr04.root
        tag="--tag Data13TeV_JetHT_2017D" #MC13TeV_AJJ_EWK_INT_LO_mjj500_dr04"

        #input=${eosdir}/Data13TeV_SinglePhoton_2017F/MergedMiniEvents_0_ext0.root
        #output=Data13TeV_SinglePhoton_2017F.root
        #tag="--tag Data13TeV_SinglePhoton_2017F"

        #input=${eosdir2018}/Data13TeV_EGamma_2018A/MergedMiniEvents_0_ext0.root
        #output=Data13TeV_EGamma_2018A.root
        #tag="--tag Data13TeV_EGamma_2018A"

	python scripts/runLocalAnalysis.py \
            -i ${input} -o ${output} ${tag} \
            --njobs 1 -q local --genWeights genweights_${githash}.root \
            --era era2017 --debug -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts --CR;

        #--debug --mvatree \
        ;;

    SEL )
	extraOpts="" #"--mvatree"
	python scripts/runLocalAnalysis.py \
	    -i ${eosdir} \
            -o ${outdir}/${githash}/${EXTRA} \
            --farmappendix ${githash} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --only Data13TeV --runSysts ${extraOpts};
	;;
        #extraOpts=" --mvatree"

    SEL2018 )
	python scripts/runLocalAnalysis.py -i ${eosdir2018} \
            -o ${outdir}/${githash2018}/raw2018 \
            -q ${queue} \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts;
	;;

    SELJETHT )
	extraOpts=" --CR"
	python scripts/runLocalAnalysis.py \
	    -i ${eosdirJetHT} \
            -o ${outdir}/${githashJetHT}/rawJetHT \
            --farmappendix ${githashJetHT} \
            -q ${queue} --genWeights genweights_${githashJetHT}.root \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --only JetHT.json --runSysts ${extraOpts};
	;;

    SELWEIGHTED )
	python scripts/runLocalAnalysis.py -i ${eosdir} \
            -o ${outdir}/weighted --flag 1 \
            -q ${queue} \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts;
        ;;

    MERGE )
        gh=${githash}
        if [[ "${EXTRA}" = *"2018"* ]]; then
            gh=${githash2018}
        fi
	if [[ "${EXTRA}" = *"JetHT"* ]]; then
	    gh=${githashJetHT}
	fi
	./scripts/mergeOutputs.py ${outdir}/${gh}/${EXTRA};
	;;

    PLOT )
        json=data/era2017/vbf_samples.json;
        lumi=${vbflumi}        
        gh=${githash}/
        if [[ "${EXTRA}" = *"2018"* ]]; then
            json=data/era2018/vbf_samples.json;
            lumi=${fulllumi2018}
            vbflumi=${lumi}
            gh=${githash2018}
        fi
        kFactors="--procSF MC13TeV_QCDEM_15to20:1.26,MC13TeV_QCDEM_20to30:1.26,MC13TeV_QCDEM_30to50:1.26,MC13TeV_QCDEM_50to80:1.26,MC13TeV_QCDEM_80to120:1.26,MC13TeV_QCDEM_120to170:1.26,MC13TeV_QCDEM_170to300:1.26,MC13TeV_QCDEM_300toInf:1.26,MC13TeV_GJets_HT40to100:1.26,MC13TeV_GJets_HT100to200:1.26,MC13TeV_GJets_HT200to400:1.26,MC13TeV_GJets_HT600toInf:1.26"
	commonOpts="-i ${outdir}/${gh}/${EXTRA} --puNormSF puwgtctr -l ${lumi}  --saveLog --mcUnc ${lumiUnc} --lumiSpecs VBFA:${vbflumi},OfflineVBFA:${fulllumi}"
	python scripts/plotter.py ${commonOpts} -j ${json} ${kFactors} --only HighMJJ -O ${outdir}/${githash}/${EXTRA}/plots/
	#python scripts/plotter.py ${commonOpts} -j ${json} ${kFactors} --only evcount --saveTeX --o ${outdir}/${githash}/${EXTRA}/plots/evcount_plotter.root;
        if [[ "${EXTRA}" != *"2018"* ]]; then
            python scripts/plotter.py ${commonOpts}  -j data/era2017/vbf_signal_samples.json --only HighPtA_ -O ${outdir}/${githash}/${EXTRA}/plots_signal/ --noStack;
            continue
            #python scripts/plotter.py ${commonOpts}  -j data/era2017/gjets_samples.json --only HighPtA_mjj -O ${outdir}/${githash}/${EXTRA}/plots_gjets/ --noStack;
            trigOpts="-i ${outdir}/${gh}/${EXTRA} --puNormSF puwgtctr -l ${vbflumi}  --saveLog --mcUnc ${lumiUnc} ${kFactors} --only HighPtOfflineVBFA,HighPtVBFA"
            #python scripts/plotter.py ${trigOpts} -j data/era2017/vbf_samples_2017F.json -O ${outdir}/${githash}/${EXTRA}/plots_trigger;
        fi
        ;;
       
    RATIO )
        #python test/analysis/computeVBFRatios.py -t \
        #    -i ${outdir}/${githash}/${EXTRA}/plots_trigger/plotter.root,${outdir}/${githash2018}/raw2018/plots/plotter.root \
        #    --titles "2017","2018" \
        #    -o ${outdir}/${githash}/${EXTRA}/plots_trigger/trigger_ratio_plotter.root

        python test/analysis/computeVBFRatios.py \
            -i ${outdir}/${githash}/${EXTRA}/plots/plotter.root,${outdir}/${githash2018}/raw2018/plots/plotter.root \
            --titles "2017","2018" \
            -o ${outdir}/${githash}/${EXTRA}/plots/ratio_plotter.root
        
	;;

    WWW )

        plotList=(plots_signal) # plots plots_trigger plots_signal plots_gjets)
        gh=${githash}
        if [[ "${EXTRA}" = *"2018"* ]]; then
            gh=${githash2018}
        fi

        for p in ${plotList[@]}; do
            pdir=${outdir}/${gh}/${EXTRA}/${p}
            if [ -d ${pdir} ]; then
                fdir=${wwwdir}/${gh}/${EXTRA}/${p}
	        mkdir -p ${fdir}
	        cp ${pdir}/*.{png,pdf,dat} ${fdir};
	        cp test/index.php ${fdir};
                echo "Check plots in ${fdir} for ${p}"
            fi
        done
	;;
esac
