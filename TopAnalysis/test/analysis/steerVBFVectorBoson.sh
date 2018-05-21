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
eosdir2018=/store/cmst3/group/top/RunIISpring18/4ad3a45
fulllumi=41367
vbflumi=7661
fulllumi2018=2300
lumiUnc=0.025
outdir=${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/analysis/VBFVectorBoson
wwwdir=~/www/VBFVectorBoson


RED='\e[31m'
NC='\e[0m'
case $WHAT in

    TESTSEL )
        input=${eosdir}/MC13TeV_DY50toInf/MergedMiniEvents_4_ext0.root
        output=MC13TeV_DY4Jets50toInf.root
        tag="--tag MC13TeV_DY50toInf"

        input=${eosdir2018}/Data13TeV_EGamma_2018A/MergedMiniEvents_0_ext0.root
        output=Data13TeV_EGamma_2018A.root
        tag="--tag Data13TeV_EGamma_2018A"

	python scripts/runLocalAnalysis.py \
            -i ${input} -o ${output} ${tag} \
            --njobs 1 -q local --debug --mvatree \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts;
        ;;

    SEL )
	python scripts/runLocalAnalysis.py -i ${eosdir} \
            -o ${outdir}/raw \
            --only data/era2017/vbf_samples.json \
            -q ${queue} \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts;
         # --mvatree \
	;;

    SEL2018 )
	python scripts/runLocalAnalysis.py -i ${eosdir2018} \
            -o ${outdir}/raw2018 \
            -q ${queue} \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts;
	;;

    SELWEIGHTED )
	python scripts/runLocalAnalysis.py -i ${eosdir} \
            -o ${outdir}/weighted --flag 1 \
            -q ${queue} \
            --era era2017 -m VBFVectorBoson::RunVBFVectorBoson --ch 0 --runSysts;
        ;;

    MERGE )
	./scripts/mergeOutputs.py ${outdir}/${EXTRA};
	;;

    PLOT )
#	commonOpts="-i ${outdir} --puNormSF puwgtctr -l 1  --saveLog"
#	python scripts/plotter.py ${commonOpts} -j data/era2017/vbf_samples.json  --noStack --skip TT,ZZ,WW,WZ,Single,QCD,GJets,DY1Jets,DY2Jets,DY3Jets,DY4Jets,DY50toInf_HT -O ${outdir}/plots_DY ;

        json=data/era2017/vbf_samples.json;
        lumi=${fulllumi}        
        if [[ "${EXTRA}" = *"2018"* ]]; then
            json=data/era2018/vbf_samples.json;
            lumi=${fulllumi2018}
            vbflumi=${lumi}
        fi

	commonOpts="-i ${outdir}/${EXTRA} --puNormSF puwgtctr -l ${fulllumi}  --saveLog --mcUnc ${lumiUnc} --lumiSpecs VBFA:${vbflumi},OfflineVBFA:${fulllumi}"
	python scripts/plotter.py ${commonOpts} -j ${json};
        ;;
    RATIO )
        python test/analysis/computeVBFRatios.py -t \
            -i ${outdir}/raw/plots/plotter.root,${outdir}/raw2018/plots/plotter.root \
            --titles "2017","2018" \
            -o ${outdir}/raw/plots/trigger_ratio_plotter.root
	;;

    WWW )
        pdir=${outdir}/${EXTRA}/plots
        fdir=${wwwdir}/${EXTRA}
	mkdir -p ${fdir}
	cp ${pdir}/*.{png,pdf} ${fdir};
	cp test/index.php ${fdir};
        echo "Check plots in ${fdir}"
	;;
esac
