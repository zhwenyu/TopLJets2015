import os

def getEraConfiguration(era,isData):
    globalTags = {
        'era2016':('80X_mcRun2_asymptotic_2016_TrancheIV_v7', '80X_dataRun2_2016SeptRepro_v6'),
        'era2017':('94X_mc2017_realistic_v10',                '94X_dataRun2_ReReco_EOY17_v2'),
        'era2018':('101X_dataRun2_Prompt_v9',                 '101X_dataRun2_Prompt_v9'),
        }
    jecTags    = {
        'era2016':('Summer16_23Sep2016V4',  'Summer16_23Sep2016AllV4'),
        'era2017':('Fall17_17Nov2017_V6',   'Fall17_17Nov2017BCDEF_V6'),
        'era2018':('Fall17_17Nov2017_V6',   'Fall17_17Nov2017BCDEF_V6'),
        }
    egmData={'era2017':['RecoEgamma/ElectronIdentification/data/Fall17',
                        'RecoEgamma/PhotonIdentification/data/Fall17/',
                        'EgammaAnalysis/ElectronTools/data/ScalesSmearings/']}

    idx=1 if isData else 0
    globalTag = globalTags[era][idx]
    print '[EraConfig] with globalTag=',globalTag

    jecTag,jecDB=None,None
    if era in jecTags:
        jecTag = jecTags[era][idx]
        jecTag_pf = 'DATA' if isData else 'MC'
        jecDB  = 'jec_DATA.db'  if isData else 'jec_MC.db'
        if not os.path.isfile(jecDB) and not os.path.islink(jecDB):
            print 'Creating symbolic link to DB correction DB'
            os.system('ln -s ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/%s/%s_%s.db %s'%(era,jecTag,jecTag_pf,jecDB))
    else:
        print '[EraConfig] No JEC tags set for',era

    if era in egmData:
        for d in egmData[era]:
            symbLinkDir='%s/src/%s'%(os.environ['CMSSW_BASE'],d)
            os.system('mkdir -p %s'%symbLinkDir) 
            ext_d='%s/external/%s/data/%s'%(os.environ['CMSSW_BASE'],os.environ['SCRAM_ARCH'],d)
            for f in os.listdir(ext_d):
                full_f=os.path.join(ext_d,f)
                if not os.path.isfile(full_f): continue            
                if os.path.isfile('%s/%s'%(symbLinkDir,f)) : continue
                os.system('cp -v %s %s/%s'%(full_f,symbLinkDir,f))
    else:
        print '[EraConfig] No EGM data specs set for',era

    return globalTag, jecTag, jecDB
    
