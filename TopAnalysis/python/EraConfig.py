import os

def getEraConfiguration(era,isData):
    globalTags = {'era2016':('80X_mcRun2_asymptotic_2016_miniAODv2_v1', '80X_dataRun2_2016SeptRepro_v5')  }
    jecTags    = {'era2016':('Spring16_23Sep2016V2',                    'Spring16_23Sep2016AllV2')}
    muonCorrs  = {'era2016':('RoccoR_13tev',                            'RoccoR_13tev')}

    idx=1 if isData else 0
    globalTag = globalTags[era][idx]
    jecTag = jecTags[era][idx]
    jecTag_pf = 'DATA' if isData else 'MC'
    jecDB  = 'jec_DATA.db'  if isData else 'jec_MC.db'
    muonDB = 'muon_DATA.txt' if isData else 'muon_MC.txt'
    if not os.path.isfile(jecDB) and not os.path.islink(jecDB):
        print 'Creating symbolic link to DB correction DB'
        os.system('ln -s ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s_%s.db %s'%(era,jecTag,jecTag_pf,jecDB))
    if not os.path.isfile(muonDB) and not os.path.islink(muonDB):
        print 'Creating symbolic link to muon correction txt file'
        os.system('ln -s ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s.txt %s'%(era,muonCorrs[era][isData],muonDB))

    return globalTag, jecTag, jecDB, muonDB
    
