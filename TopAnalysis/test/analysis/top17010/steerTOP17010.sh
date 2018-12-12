#!/bin/bash

ERA=2016

while getopts "o:e:q:y:" opt; do
    case "$opt" in
        o) WHAT=$OPTARG
            ;;
        y) ERA=$OPTARG
            ;;
    esac
done

if [ -z "$WHAT" ]; then 
    echo "steerTopWidth.sh -o <SEL/MERGE/...> [-y 2016/7]";
    echo "        RESOL        - preliminary resolution study"
    echo "        TESTSEL      - test selection locally"
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGE        - merge output (if given \"extra\" is appended to the directory)"
    echo "        PLOT         - make plots (if given \"extra\" is appended to the directory)"
    echo "        WWW          - move plots to web-based (if given \"extra\" is appended to the directory)"
    exit 1; 
fi

githash=0c522df
eosdir=/store/cmst3/group/top/RunIIReReco/2016/${githash}
fulllumi=35505
lumiUnc=0.025
if [[ ${ERA} = "2017" ]]; then
    githash=3129835
    eosdir=/store/cmst3/group/top/RunIIReReco/${githash}
    fulllumi=41367
fi
echo "Selection adapted to YEAR=${ERA}, inputs from ${eosdir}"

queue=workday
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/top17010
json=test/analysis/top17010/samples_${ERA}.json
systJson=test/analysis/top17010/syst_samples_${ERA}.json
wwwdir=${HOME}/www/top17010

mkdir -p ${outdir}

RED='\e[31m'
NC='\e[0m'
case $WHAT in

    RESOL )
        tag=MC13TeV_${ERA}_TTJets_m1755
        python test/analysis/top17010/getMlbResolution.py /eos/cms/${eosdir}/${tag}/Chunk_0_ext0.root
        mkdir -p ${outdir}/plots
        mv mlbresol_*.{png,pdf} ${outdir}/plots
        ;;

    TESTSEL )               
        tag=MC13TeV_${ERA}_TTJets
        input=${eosdir}/${tag}/Chunk_0_ext0.root
        output=${tag}.root 
	python scripts/runLocalAnalysis.py \
            -i ${input} -o ${output} --tag ${tag} --only ${json} --flag 0\
            --njobs 1 -q local --genWeights genweights_${githash}.root \
            --era era${ERA} -m TOP17010::TOP17010 --ch 0 --runSysts; --debug;
        ;;

    SEL )
        
	python scripts/runLocalAnalysis.py \
	    -i ${eosdir} --only ${json}\
            -o ${outdir}/${githash}/${EXTRA} \
            --farmappendix ${githash} \
            -q ${queue} --genWeights genweights_${githash}.root \
            --era era${ERA} -m  TOP17010::TOP17010 --ch 0 --runSysts;
	;;

    MERGE )
	./scripts/mergeOutputs.py ${outdir}/${githash}/${EXTRA};
	;;

    PLOT )
	#FIXME
        json=data/era${ERA}/vbf_samples.json;
	syst_json=data/${ERA}/vbf_syst_samples.json;
	plotOutDir=${outdir}/${githash}/${EXTRA}/plots/
        kFactors="--procSF MC13TeV_era${ERA}_QCDEM_15to20:1.26,MC13TeV_era${ERA}_QCDEM_20to30:1.26,MC13TeV_era${ERA}_QCDEM_30to50:1.26,MC13TeV_era${ERA}_QCDEM_50to80:1.26,MC13TeV_era${ERA}_QCDEM_80to120:1.26,MC13TeV_era${ERA}_QCDEM_120to170:1.26,MC13TeV_era${ERA}_QCDEM_170to300:1.26,MC13TeV_era${ERA}_QCDEM_300toInf:1.26,MC13TeV_era${ERA}_GJets_HT40to100:1.26,MC13TeV_era${ERA}_GJets_HT100to200:1.26,MC13TeV_era${ERA}_GJets_HT200to400:1.26,MC13TeV_era${ERA}_GJets_HT600toInf:1.26"
	commonOpts="-i ${outdir}/${githash}/${EXTRA} --puNormSF puwgtctr -l ${fulllumi} --saveLog --mcUnc ${lumiUnc} --lumiSpecs LowVPtLowMJJA:${vbflumi},LowVPtHighMJJA:${vbflumi}"
	python scripts/plotter.py ${commonOpts} -j ${json} --only HighMJJ,LowMJJ ${kFactors}
        python scripts/plotter.py ${commonOpts} -j ${json} --only evcount ${kFactors} --saveTeX -o evcout_plotter.root
	python scripts/plotter.py ${commonOpts} -j ${syst_json} ${kFactors} --only HighMJJ,LowMJJ --silent -o syst_plotter.root
        ;;
    
    WWW )
        #FIXME
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
