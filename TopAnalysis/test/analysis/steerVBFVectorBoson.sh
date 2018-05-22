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
fulllumi=41367
vbflumi=7661
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
        #input=${eosdir}/Data13TeV_SingleMuon_2017D/MergedMiniEvents_0_ext0.root
        #output=Data13TeV_SingleMuon_2017D.root
        #input=${eosdir}/MC13TeV_GJets_HT100to200_DR04/MergedMiniEvents_0_ext0.root
        #output=MC13TeV_GJets_HT100to200_DR04.root \
        #input=${eosdir}/Data13TeV_SingleMuon_2017C/MergedMiniEvents_0_ext0.root
        #output=Data13TeV_SingleMuon_2017C.root
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
	commonOpts="-i ${outdir}/${EXTRA} --puNormSF puwgtctr -l ${fulllumi}  --saveLog --mcUnc ${lumiUnc} --lumiSpecs VBFA:${vbflumi}"
	python scripts/plotter.py ${commonOpts} -j data/era2017/vbf_samples.json; 
        python test/analysis/computeVBFRatios.py -t -i ${outdir}/${EXTRA}/plots/plotter.root -o ${outdir}/${EXTRA}/plots/trigger_ratio_plotter.root

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
