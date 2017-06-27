#!/usr/bin/env python2.7

from collections import OrderedDict

def defineDataset(data):

    baseDir='root://eoscms//eos/cms/store/cmst3/group/top/HIN-17-002/'

    if 'Data8.16TeV_pPb' in data:
        dataset=OrderedDict()
        dataset['l+jets']         = (['%s/Data_PbP_EJets_GEDElectrons_WithHFAndSkimInfo_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root'%baseDir,
                                      '%s/Data_PPb_EJets_GEDElectrons_WithHFAndSkimInfo_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root'%baseDir,
                                      '%s/Data_PbP_MuJets_WithHFAndSkimInfo.root',
                                      '%s/Data_PPb_MuJets_WithHFAndSkimInfo.root'],
                                     '#4575b4', [11,13],   False)
        dataset['non-iso l+jets'] = (['%s/DataNonIso_PbP_EJets_GEDElectrons_WithHFAndSkimInfo_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root'%baseDir,
                                      '%s/DataNonIso_PPb_EJets_GEDElectrons_WithHFAndSkimInfo_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root'%baseDir,
                                      '%s/DataNonIso_PPb_MuJets_WithHFAndSkimInfo.root'%baseDir,
                                      '%s/DataNonIso_PbP_MuJets_WithHFAndSkimInfo.root'%baseDir],
                                      '#91bfdb', [1100,1300], False)
    if 'Data8.16TeV_pPb_nonsubiso' in data:
        dataset=OrderedDict()
        dataset['l+jets']         = (['%s/Data_PbP_MuJets_WithHFAndSkimInfo_NonSubractedIso.root'%baseDir,
                                      '%s/Data_PPb_MuJets_WithHFAndSkimInfo_NonSubractedIso.root'%baseDir,
                                      '%s/Data_PbP_EJets_GEDElectrons_WithHFAndSkimInfo_NonSubractedIso_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root'%baseDir,
                                      '%s/Data_PPb_EJets_GEDElectrons_WithHFAndSkimInfo_NonSubractedIso_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root'%baseDir],
                                     '#4575b4', [11,13],   False)
        dataset['non-iso l+jets'] = (['%s/DataNonIso_PPb_MuJets_WithHFAndSkimInfo_NonSubractedIso.root',
                                      '%s/DataNonIso_PbP_MuJets_WithHFAndSkimInfo_NonSubractedIso.root',
                                      '%s/DataNonIso_PbP_EJets_GEDElectrons_WithHFAndSkimInfo_NonSubractedIso_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root'%baseDir,
                                      '%s/DataNonIso_PPb_EJets_GEDElectrons_WithHFAndSkimInfo_NonSubractedIso_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root'%baseDir],
                                     '#91bfdb', [1100,1300], False)
    if 'MC8.16TeV_TTbar_pPb' in data:
        dataset=OrderedDict()

        dataset['l+jets']         = (['%s/PYQUEN_TTbar_PPb-EmbEPOS_EJets.root'%baseDir,
                                      '%s/PYQUEN_TTbar_PPb-EmbEPOS_MuJets.root'%baseDir],
                                      '#4575b4', [11,13],   True)
    if 'MC8.16TeV_TTbar_pPb_Pohweg' in data:
        dataset=OrderedDict()
        dataset['l+jets']         = (['%s/Powheg_TTbar_PPb-EmbEPOS_MuJets.root'%baseDir,
                                      '%s/Powheg_TTbar_PbP-EmbEPOS_MuJets.root'%baseDir,
                                      '%s/Powheg_TTbar_PbP-EmbEPOS_EJets.root'%baseDir,
                                      '%s/Powheg_TTbar_PPb-EmbEPOS_EJets.root'%baseDir],
                                      '#4575b4', [11,13],   True)

    if 'MC8.16TeV_WJets_pPb' in data:
        dataset=OrderedDict()
        dataset['l+jets']         = (['%s/WToMuNu_anstahll-Embedded_PileUp_POWHEG_WToMuNu_Plus_CT14_EPPS16_8160GeV_pythia8_reverse_RECO.root'%baseDir,
                                      '%s/WToMuNu_anstahll-Embedded_PileUp_POWHEG_WToMuNu_Minus_CT14_EPPS16_8160GeV_pythia8_reverse_RECO.root'%baseDir,
                                      '%s/WToMuNu_anstahll-Embedded_PileUp_POWHEG_WToMuNu_Plus_CT14_EPPS16_8160GeV_pythia8_RECO.root'%baseDir,
                                      '%s/WToMuNu_anstahll-Embedded_PileUp_POWHEG_WToMuNu_Minus_CT14_EPPS16_8160GeV_pythia8_RECO.root'%baseDir],
                                     '#4575b4', [11,13],   False)

    if 'MC8TeV_WJets_pp' in data:
        dataset=OrderedDict()
        dataset['l+jets'] = ( ['%s/MC8TeV_W3Jets.root'%baseDir],
                               '#4575b4',
                               [11,13],#[11,13],
                               False )
    if data=='Data8TeV_pp':
        dataset=OrderedDict()
        dataset['l+jets'] = ( ['%s/Data8TeV_SingleMu2012A.root'%baseDir
                                #'%s/Data8TeV_SingleMu2012B.root'%baseDir,
                                #'%s/Data8TeV_SingleMu2012C.root'%baseDir,
                                #'%s/Data8TeV_SingleMu2012D.root%baseDir'
                                ],
                                 '#4575b4',
                                 [11,13],
                                 False )

    if data=='MC8TeV_TTJets':
        dataset=OrderedDict()
        dataset['l+jets'] = ( ['%s/MC8TeV_TTJets_MSDecays_172v5.root'%baseDir],
                              '#4575b4',
                              [13], #[11,13]
                              True )

    return dataset
