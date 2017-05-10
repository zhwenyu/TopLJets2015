#!/bin/bash

WHAT=$1; 

case $WHAT in

    TESTSEL)

	for sample in MC8.16TeV_WJets_pPb MC8TeV_WJets_pp MC8.16TeV_TTbar_pPb; do
	    for wjjOrder in drjj mjj sumpt; do
		for thadOrder in dr dm2tlep dm2pdg; do
	            python2.7 prepareWorkspace.py  -d ${sample} --wjjOrder ${wjjOrder} --thadOrder ${thadOrder};
	            python2.7 checkMatchingEfficiency.py plots/${sample}/controlplots.root > plots/${sample}/algoeff_${wjjOrder}_${thadOrder}.tex
	            mv plots/${sample}/controlplots.root plots/${sample}/controlplots_${wjjOrder}_${thadOrder}.root;
		done
	    done

	    rm plots/${sample}/*.{png,pdf}
	    
	    python2.7 checkMatchingDistributions.py "min#DeltaR(jj)":plots/${sample}/controlplots_drjj_dm2tlep.root \
	        "max#Sigmap_{T}(j)":plots/${sample}/controlplots_sumpt_dm2tlep.root \
	        "min|M(jj)-M(W)|":plots/${sample}/controlplots_mjj_dm2tlep.root
	    mkdir plots/${sample}/w
	    mv *.{png,pdf} plots/${sample}/w
	    
	    
	    python2.7 checkMatchingDistributions.py "min|m_{thad}-m_{tlep}|":plots/${sample}/controlplots_drjj_dm2tlep.root \
		"min|m_{thad}_m_{t}|":plots/${sample}/controlplots_sumpt_dm2pdg.root \
		"min#DeltaR(b,W)":plots/${sample}/controlplots_mjj_dr.root
	    mkdir plots/${sample}/top
	    mv *.{png,pdf} plots/${sample}/top
	done
	;;
    
    SEL)
	common="--wjjOrder drjj --thadOrder dm2tlep"

	#JER/JES variations
	sample=MC8.16TeV_TTbar_pPb
	for jer in -1 1; do
	    python2.7 prepareWorkspace.py  -d ${sample} ${common} --smearJER ${jer};
	    mv workspace_MC8.16TeV_TTbar_pPb.root workspace_MC8.16TeV_TTbar_pPb_jer${jer}.root
	done
	for jes in 0.964 1.036; do
	    python2.7 prepareWorkspace.py  -d ${sample} ${common} --shiftJES ${jes};
	    mv workspace_MC8.16TeV_TTbar_pPb.root workspace_MC8.16TeV_TTbar_pPb_jes${jes}.root
	done
	#nominal
	for sample in MC8.16TeV_TTbar_pPb MC8.16TeV_WJets_pPb; do
	    python2.7 prepareWorkspace.py  -d ${sample} ${common};
	done
	;;

    SELDATA )
	#nominal
	for sample in Data8.16TeV_pPb_mu Data8.16TeV_pPb_e; do
	    python2.7 prepareWorkspace.py  -d ${sample} ${common};
	done
	;;

    PARAM)

	a=(`ls workspace_MC8.16TeV_TTbar_pPb*root`)
	for i in ${a[@]}; do
	    python2.7 parameterizeMCShapes.py ${i};	    
	    outDir=plots/MC8.16TeV_TTbar_pPb/${i/.root/}	    
	    mkdir ${outDir}
	    b=(`ls *.{png,pdf}`)
	    for i in ${b[@]}; do
		f=${i/./_ttbar.}
		mv -v ${i} ${outDir}/${f}; 
	    done
	done

	python2.7 parameterizeMCShapes.py workspace_MC8.16TeV_WJets_pPb.root pdf_workspace_MC8.16TeV_TTbar_pPb.root;
	a=(`ls *.{png,pdf}`)
	for i in ${a[@]}; do
	    f=${i/./_w.}
	    mv -v ${i} plots/MC8.16TeV_WJets_pPb/${f}; 
	done
	;;

    FITS)
	for onlyResonant in 0 1 2; do
	    for l in e mu; do
		commonOpts="--onlyResonant ${onlyResonant} -s pdf_workspace_MC8.16TeV_TTbar_pPb.root"
		python2.7 runDataFit.py ${commonOpts} -i workspace_Data8.16TeV_pPb_${l}.root -o plots/Data8.16TeV_pPb_${l}/;
	    done
	done

	;;

    SYSTS)
	for l in e mu; do
	    for syst in jer-1 jer1 jes0.964 jes1.036; do
		commonOpts="--onlyResonant 0 -s pdf_workspace_MC8.16TeV_TTbar_pPb_${syst}.root";
		mkdir -p plots/Data8.16TeV_pPb_${l}_${syst}/
		python2.7 runDataFit.py ${commonOpts} -i workspace_Data8.16TeV_pPb_${l}.root -o plots/Data8.16TeV_pPb_${l}_${syst}/;
	    done
	done
	;;

    FINALCOMP )
	common="--wjjOrder drjj --thadOrder dm2tlep"
	for l in e mu; do
	    python2.7 prepareWorkspace.py  -d MC8.16TeV_TTbar_pPb_${l} ${common};
	done
	python2.7 combineFinalControlPlots.py
	;;
esac
