#!/bin/bash

while getopts "o:e:q:y:d:" opt; do
    case "$opt" in
        o) WHAT=$OPTARG
            ;;
        e) EXTRA=$OPTARG
            ;;
        q) QCD=$OPTARG
            ;;
        y) ERA=$OPTARG
            ;;
    esac
done

if [ -z "$WHAT" ]; then 
    echo "steerVBFVectorBoson.sh -o <SEL/MERGE/...> [-e extra -q QCD -y 2016/7]";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGE        - merge output (if given \"extra\" is appended to the directory)"
    echo "        PLOT         - make plots (if given \"extra\" is appended to the directory)"
    echo "        NLOTFACTORS  - compute NLO/LO in DY and test how it looks after applying to gamma+jets"
    echo "        TRIGEFF      - trigger efficiency"
    echo "        BDTTRANSFORM - derive BDT transformation for uniform background yields"
    echo "        WWW          - move plots to web-based (if given \"extra\" is appended to the directory)"
    exit 1; 
fi

githash=3129835
eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
fulllumi=41367
vbflumi=7661
lumiUnc=0.025
if [[ ${ERA} == "2016" ]]; then
    githash=0c522df
    eosdir=/store/cmst3/group/top/RunIIReReco/2016/${githash}
    fulllumi=35900
    vbflumi=28000
fi

echo "Selection adapted to YEAR=${ERA}"

#to run locally use local as queue + can add "--njobs 8" to use 8 parallel jobs
queue=workday
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBFVectorBoson
wwwdir=~/www/VBFVectorBoson


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
               
        json=data/era${ERA}/vbf_samples.json
        #tag=Data13TeV_2017C_SingleMuon
        #if [[ ${ERA} == "2016" ]]; then
        #    tag=MC13TeV_2016_TTJets
        #fi
        tag=MC13TeV_${ERA}_EWKAJJ
        input=${eosdir}/${tag}/Chunk_1_ext0.root
        output=${tag}.root 

	python scripts/runLocalAnalysis.py \
            -i ${input} -o ${output} --tag ${tag} --only ${json} --mvatree\
            --njobs 1 -q local --genWeights genweights_${githash}.root \
            --era era${ERA} -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts --debug;

        #--debug --mvatree \
        ;;


    TESTSELTRIGEFF )
               
        json=data/era${ERA}/vbf_samples.json
        tag=Data13TeV_${ERA}C_SingleMuon
        input=${eosdir}/${tag}/Chunk_1_ext0.root
        output=${tag}.root 

	python scripts/runLocalAnalysis.py \
            -i ${input} -o ${output} --tag ${tag} --only ${json} --mvatree\
            --njobs 1 -q local --genWeights genweights_${githash}.root \
            --era era${ERA} -m PhotonTrigEff::RunPhotonTrigEff --ch 0 --runSysts --debug;

        ;;

    SEL )
	##### NOTE: There are three options here:
        ### --mvatree: to store trees for BDT training in signal region
        ### --CR     : gives a control region to evaluate fake rates in the photon data samples
        ### --SRfake : gives the distributions of fakes, normalised based on fake rates

        json=data/era${ERA}/vbf_samples.json,data/era${ERA}/vbf_syst_samples.json
	if [[ -z ${EXTRA} ]]; then
	    echo "Making trees ... "
	    extraOpts=" --mvatree"
	    json=data/era${ERA}/vbf_trees.json
	    EXTRA="MVATrees"
        fi
	python scripts/runLocalAnalysis.py \
	    -i ${eosdir} --only ${json} \
            -o ${outdir}/${githash}/${EXTRA} \
            --farmappendix ${githash} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era${ERA} -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts ${extraOpts};
	;;

    SELTRIGEFF )
        json=data/era${ERA}/vbf_samples.json 
	python scripts/runLocalAnalysis.py \
	    -i ${eosdir} --only ${json}\
            -o ${outdir}/trig/${githash}/${EXTRA} \
            --farmappendix trig${githash} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era${ERA} -m PhotonTrigEff::RunPhotonTrigEff --ch 0 --runSysts ${extraOpts};
	;;


    SELJETHT )
	json=data/era${ERA}/JetHT.json;
	extraOpts=" --CR"
	echo ${QCD} 
	if [[ ${QCD} == "QCDTemp" ]]; then
	    echo 'I do QCD Template photon selection'
	    extraOpts=${extraOpts}" --QCDTemp"
	fi
        if [[ ${QCD} == "SRfake" ]]; then
            echo 'I do SRfake photon selection'
            extraOpts=" --SRfake"
        fi
	python scripts/runLocalAnalysis.py \
	    -i ${eosdir} --only ${json} \
            -o ${outdir}/${githash}/${EXTRA}${QCD} \
            --farmappendix ${githash}${EXTRA}${QCD} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era${ERA} -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts ${extraOpts} --skipexisting;
	;;


    MERGE )
	./scripts/mergeOutputs.py ${outdir}/${githash}/${EXTRA};
	;;

    MERGETRIGEFF )
	./scripts/mergeOutputs.py ${outdir}/trig/${githash}/${EXTRA};
	;;

    PLOT )
	
        json=data/era${ERA}/vbf_samples.json;
	syst_json=data/era${ERA}/vbf_syst_samples.json;
        gjets_json=data/era${ERA}/gjets_samples.json;
	plotOutDir=${outdir}/${githash}/${EXTRA}/plots/
        kFactors="--procSF MC13TeV_era${ERA}_QCDEM_15to20:1.26,MC13TeV_era${ERA}_QCDEM_20to30:1.26,MC13TeV_era${ERA}_QCDEM_30to50:1.26,MC13TeV_era${ERA}_QCDEM_50to80:1.26,MC13TeV_era${ERA}_QCDEM_80to120:1.26,MC13TeV_era${ERA}_QCDEM_120to170:1.26,MC13TeV_era${ERA}_QCDEM_170to300:1.26,MC13TeV_era${ERA}_QCDEM_300toInf:1.26,MC13TeV_era${ERA}_GJets_HT40to100:1.26,MC13TeV_era${ERA}_GJets_HT100to200:1.26,MC13TeV_era${ERA}_GJets_HT200to400:1.26,MC13TeV_era${ERA}_GJets_HT400to600:1.26,MC13TeV_era${ERA}_GJets_HT600toInf:1.26"
	commonOpts="-i ${outdir}/${githash}/${EXTRA} --puNormSF puwgtctr -l ${fulllumi} --saveLog --mcUnc ${lumiUnc} --lumiSpecs LowVPtLowMJJA:${vbflumi},LowVPtHighMJJA:${vbflumi}"
        #python scripts/plotter.py ${commonOpts} -j ${gjets_json} --noStack --only A_
	python scripts/plotter.py ${commonOpts} -j ${json} --only HighMJJ,LowMJJ ${kFactors}
        #python scripts/plotter.py ${commonOpts} -j ${json} --only evcount ${kFactors} --saveTeX -o evcout_plotter.root
	#python scripts/plotter.py ${commonOpts} -j ${syst_json} ${kFactors} --only HighMJJ,LowMJJ --silent -o syst_plotter.root
        ;;
    
    NLOTFACTORS )
        #transfer factors
        plotOutDir=${outdir}/${githash}/${EXTRA}/plots/
        commonOpts="-p ${plotOutDir}/plotter.root -s ${plotOutDir}/syst_plotter.root -o ${plotOutDir}"
        python test/analysis/computeTransferFactor.py ${commonOpts} --var leadpt     --binList 50,100,150,200,250,300,400,500
        python test/analysis/computeTransferFactor.py ${commonOpts} --var centraleta --binList 0,0.4,0.8,1.2,2,3,5
        python test/analysis/computeTransferFactor.py ${commonOpts} --var forwardeta --binList 0,0.8,1.2,2,3,5
        python test/analysis/computeTransferFactor.py ${commonOpts} --var detajj     --binList 0,0.8,1.6,2.4,3.2,4,5,8
        python test/analysis/computeTransferFactor.py ${commonOpts} --var dphijj     --rebin 2
        python test/analysis/computeTransferFactor.py ${commonOpts} --var mjj        --binList 500,600,700,800,900,1000,1250,1500,2000,3000,4000
        python test/analysis/computeTransferFactor.py ${commonOpts} --var vpt        --binList 75,100,150,200,250,300,350,400,550
        python test/analysis/computeTransferFactor.py ${commonOpts} --var vbfmva     --rebin 5
        ;;

    TRIGEFF )
        #trigger efficiencies
	inDir=${outdir}/trig/${githash}/${EXTRA}
        python test/analysis/computeTriggerEff.py ${inDir} ${ERA};
        #python test/analysis/computeVBFTriggerEff.py -p ${plotOutDir}/plotter.root -o ${plotOutDir};
        ;;

    BDTTRANSFORM )
        #python test/analysis/VBF_weights/getInverseCDF.py
        python test/analysis/VBF_weights/getInverseCDFFromPlotter.py ${outdir}/${githash}/${EXTRA}/plots/plotter.root
        cp -v inverse_cdfs.root test/analysis/VBF_weights/inverse_cdfs.root
        ;;

    WWW )

        plotList=(plots)
        gh=${githash}
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
