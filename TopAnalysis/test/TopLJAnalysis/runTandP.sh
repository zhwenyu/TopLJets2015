f=(BWResCBExp) #CB2Exp)
d=(True False)
for i in ${f[@]}; do 
    for j in ${d[@]}; do                        
        continue
        cmsRun runTandP_HLT.py data=${j} func=${i} & 
	cmsRun runTandP_IDISO.py data=${j} func=${i} &  
	cmsRun runTandP_HLTIDISO.py data=${j} func=${i} & 
    done 

    python drawTandPResults.py hltfits_pPb_8.16TeV_1_BWResCBExp.root:"Data",hltfits_pPb_8.16TeV_0_BWResCBExp.root:"MC" ~/www/LJets-pPb/HLT_mu
    python drawTandPResults.py idisofits_pPb_8.16TeV_1_BWResCBExp.root:"Data",idisofits_pPb_8.16TeV_0_BWResCBExp.root:"MC" ~/www/LJets-pPb/IDISO_mu
    python drawTandPResults.py hltidisofits_pPb_8.16TeV_1_BWResCBExp.root:"Data",hltidisofits_pPb_8.16TeV_0_BWResCBExp.root:"MC" ~/www/LJets-pPb/HLTIDISO0.4_mu

done

