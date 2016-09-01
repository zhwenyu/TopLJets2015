#!/bin/bash

WHAT=$1; 
ERA=$2
if [ "$#" -ne 2 ]; then 
    echo "steerTOPWidthAnalysis.sh <SEL/MERGESEL/PLOTSEL/WWWSEL/ANA/MERGE/BKG/PLOT/WWW> <ERA>";
    echo "        SEL          - launches selection jobs to the batch, output will contain summary trees and control plots"; 
    echo "        MERGESEL     - merge the output of the jobs";
    echo "        PLOTSEL      - runs the plotter tool on the selection";
    echo "        WWWSEL       - moves the plots to an afs-web-based area";
    echo "        ANA          - analyze the selected events";
    echo "        MERGE        - merge the output of the analysis jobs";
    echo "        BKG          - estimate DY scale factor from data";
    echo "        PLOT         - runs the plotter tool on the analysis outputs";
    echo "        WWW          - moves the analysis plots to an afs-web-based area";
    echo " "
    echo "        ERA          - era2015/era2016";
    exit 1; 
fi

export LSB_JOB_REPORT_MAIL=N


queue=2nw
githash=8db9ad6
lumi=12870
lumiSpecs="--lumiSpecs EE:11391"
lumiUnc=0.062
whoami=`whoami`
myletter=${whoami:0:1}
eosdir=/store/cmst3/user/psilva/LJets2016/${githash}
summaryeosdir=/store/cmst3/group/top/summer2016/TopWidth_${ERA}_ichep
case $ERA in
    era2015)
	githash=8c1e7c9;
	lumi=2267.84
	lumiSpecs=""
	lumiUnc=0.027
	eosdir=/store/cmst3/user/psilva/LJets2015/${githash}
	summaryeosdir=/store/cmst3/group/top/summer2016/TopWidth_${ERA}_notrig 
        #/store/cmst3/group/top/summer2016/TopWidth_${ERA}
	;;
esac
COMBINERELEASE=~/scratch0/CMSSW_7_4_7/src/
outdir=/afs/cern.ch/work/${myletter}/${whoami}/TopWidth_${ERA}
wwwdir=~/www/TopWidth_${ERA}


RED='\e[31m'
NC='\e[0m'
case $WHAT in
    SEL )
	python scripts/runLocalAnalysis.py -i ${eosdir} -q ${queue} -o ${outdir} --era ${ERA} -m TOP-16-019::RunTop16019 --ch 0;
	;;
    MERGESEL )
	./scripts/mergeOutputs.py ${outdir} True;	
	;;
    PLOTSEL )
	python scripts/plotter.py -i ${outdir} --puNormSF puwgtctr  -j data/${ERA}/samples.json -l ${lumi} ${lumiSpecs} --saveLog --mcUnc ${lumiUnc};	
	;;
    WWWSEL )
	mkdir -p ${wwwdir}/sel
	cp ${outdir}/plots/*.{png,pdf} ${wwwdir}/sel
	cp test/index.php ${wwwdir}/sel
	;;
    ANA )
	python scripts/runTopWidthAnalysis.py -i ${summaryeosdir} -o ${outdir}/analysis -q ${queue};
	;;
    MERGE )
	./scripts/mergeOutputs.py ${outdir}/analysis;
	;;
    BKG )
	python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/samples.json  -l ${lumi} ${lumiSpecs} --onlyData --only mll -o dy_plotter.root;        
	python scripts/runDYRinRout.py --in ${outdir}/analysis/plots/dy_plotter.root --categs 1b,2b --out ${outdir}/analysis/plots/;
	;;
    PLOT )
	python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/samples.json      -l ${lumi} ${lumiSpecs} --mcUnc ${lumiUnc} --only count --saveTeX -o count_plotter.root --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck; 
	python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/samples.json      -l ${lumi} ${lumiSpecs} --mcUnc ${lumiUnc} --only njets,ptlb -o njets_plotter.root --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck; 
        python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/samples.json      -l ${lumi} ${lumiSpecs} --mcUnc ${lumiUnc} --onlyData --procSF DY:${outdir}/analysis/plots/.dyscalefactors.pck &
	python scripts/plotter.py -i ${outdir}/analysis  -j data/${ERA}/syst_samples.json -l ${lumi} ${lumiSpecs} --mcUnc ${lumiUnc} --silent -o syst_plotter.root;
	case $ERA in
	    era2016)
		mv ${outdir}/analysis/plots/syst_plotter.root ${outdir}/analysis/plots/syst_plotter_orig.root
		python test/TopWidthAnalysis/createtWShapeUncs.py ~/work/TopWidth_era2015/analysis/plots/plotter.root ~/work/TopWidth_era2015/analysis/plots/syst_plotter.root  ${outdir}/analysis/plots/plotter.root;
		hadd ${outdir}/analysis/plots/syst_plotter.root tW_syst_plotter.root ${outdir}/analysis/plots/syst_plotter_orig.root 
		mv tW_syst_plotter.root ${outdir}/analysis/plots/tW_syst_plotter.root;	       
	esac
	#combined plots
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb        EE1b,EE2b,MM1b,MM2b,EM1b,EM2b    ${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb        EE1b,MM1b,EM1b                   ${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py ptlb        EE2b,MM2b,EM2b                   ${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_1.0w lowptEE1b,lowptMM1b,lowptEM1b    ${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_1.0w lowptEE2b,lowptMM2b,lowptEM2b    ${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_1.0w highptEE1b,highptMM1b,highptEM1b ${outdir}/analysis/plots/plotter.root
	python test/TopWidthAnalysis/combinePlotsForAllCategories.py incmlb_1.0w highptEE2b,highptMM2b,highptEM2b ${outdir}/analysis/plots/plotter.root
        ;;
    WWW )
        mkdir -p ${wwwdir}/ana
        cp ${outdir}/analysis/plots/*.{png,pdf} ${wwwdir}/ana        
        cp test/index.php ${wwwdir}/ana
	mkdir -p ${wwwdir}/comb
        cp plots/*.{root,png,pdf} ${wwwdir}/comb 
        cp test/index.php ${wwwdir}/comb
	;;
    DATACARDS ) # to get the shapes file / datacards for the full analysis
	echo "Forcing to use CMSSW_7_4_7"
        cd ${COMBINERELEASE}
        eval `scram r -sh`
        cd -
	cardsdir=${outdir}/datacards
	mkdir -p ${cardsdir}
        rm ${cardsdir}

        for index in `seq 0 17 204` ; do
            min=$index
            max=$(($index+16))

            nohup python test/TopWidthAnalysis/createShapesFromPlotter.py \
                -s tbart,tW \
                --dists incmlb \
                -o ${outdir}/datacards/ \
                -n shapes$index \
                -i ${outdir}/analysis/plots/plotter.root \
                --systInput ${outdir}/analysis/plots/syst_plotter.root \
                --nomorph --min $min --max $max > ${cardsdir}/shapes${index}.txt & 
        done
    ;;
    
esac
