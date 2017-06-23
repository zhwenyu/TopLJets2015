#!/usr/bin/env python2.7

from collections import OrderedDict

def defineDataset(data):

    if 'Data8.16TeV_pPb' in data:
        dataset=OrderedDict()
        dataset['l+jets']         = (['/afs/cern.ch/work/g/gkrintir/public/forMarta/Data_PbP_EJets_GEDElectrons_WithHFAndSkimInfo_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/Data_PPb_EJets_GEDElectrons_WithHFAndSkimInfo_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/Data_PbP_MuJets_WithHFAndSkimInfo.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/Data_PPb_MuJets_WithHFAndSkimInfo.root'],
                                     '#4575b4', [11,13],   False)
        dataset['non-iso l+jets'] = (['/afs/cern.ch/work/g/gkrintir/public/forMarta/DataNonIso_PbP_EJets_GEDElectrons_WithHFAndSkimInfo_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/DataNonIso_PPb_EJets_GEDElectrons_WithHFAndSkimInfo_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/DataNonIso_PPb_MuJets_WithHFAndSkimInfo.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/DataNonIso_PbP_MuJets_WithHFAndSkimInfo.root'],
                                      '#91bfdb', [1100,1300], False)
    if 'Data8.16TeV_pPb_nonsubiso' in data:
        dataset=OrderedDict()
        dataset['l+jets']         = (['/afs/cern.ch/work/g/gkrintir/public/forMarta/Data_PbP_MuJets_WithHFAndSkimInfo_NonSubractedIso.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/Data_PPb_MuJets_WithHFAndSkimInfo_NonSubractedIso.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/Data_PbP_EJets_GEDElectrons_WithHFAndSkimInfo_NonSubractedIso_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/Data_PPb_EJets_GEDElectrons_WithHFAndSkimInfo_NonSubractedIso_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root'],
                                     '#4575b4', [11,13],   False)
        dataset['non-iso l+jets'] = (['/afs/cern.ch/work/g/gkrintir/public/forMarta/DataNonIso_PPb_MuJets_WithHFAndSkimInfo_NonSubractedIso.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/DataNonIso_PbP_MuJets_WithHFAndSkimInfo_NonSubractedIso.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/DataNonIso_PbP_EJets_GEDElectrons_WithHFAndSkimInfo_NonSubractedIso_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/DataNonIso_PPb_EJets_GEDElectrons_WithHFAndSkimInfo_NonSubractedIso_PAEle20_WPLoose_Gsf_or_PASinglePhoton30_Eta3p1.root'],
                                     '#91bfdb', [1100,1300], False)
    if 'MC8.16TeV_TTbar_pPb' in data:
        dataset=OrderedDict()

        dataset['l+jets']         = (['/afs/cern.ch/work/g/gkrintir/public/forMarta/PYQUEN_TTbar_PPb-EmbEPOS_EJets.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/PYQUEN_TTbar_PPb-EmbEPOS_MuJets.root'],
                                      '#4575b4', [11,13],   True)
    if 'MC8.16TeV_TTbar_pPb_Pohweg' in data:
        dataset=OrderedDict()
        dataset['l+jets']         = (['/afs/cern.ch/work/g/gkrintir/public/forMarta/Powheg_TTbar_PPb-EmbEPOS_EJets.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/Powheg_TTbar_PPb-EmbEPOS_MuJets.root'],
                                      '#4575b4', [11,13],   True)


    if 'MC8.16TeV_WJets_pPb' in data:
        dataset=OrderedDict()
        dataset['l+jets']         = (['/afs/cern.ch/work/g/gkrintir/public/forMarta/WToMuNu_anstahll-Embedded_PileUp_POWHEG_WToMuNu_Plus_CT14_EPPS16_8160GeV_pythia8_reverse_RECO.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/WToMuNu_anstahll-Embedded_PileUp_POWHEG_WToMuNu_Minus_CT14_EPPS16_8160GeV_pythia8_reverse_RECO.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/WToMuNu_anstahll-Embedded_PileUp_POWHEG_WToMuNu_Plus_CT14_EPPS16_8160GeV_pythia8_RECO.root',
                                      '/afs/cern.ch/work/g/gkrintir/public/forMarta/WToMuNu_anstahll-Embedded_PileUp_POWHEG_WToMuNu_Minus_CT14_EPPS16_8160GeV_pythia8_RECO.root'],
                                     '#4575b4', [11,13],   False)

    if 'MC8TeV_WJets_pp' in data:
        dataset=OrderedDict()
        dataset['l+jets'] = ( [#'data_pp/MC8TeV_WJets.root',
                               #'data_pp/MC8TeV_W1Jets.root',
                               #'data_pp/MC8TeV_W2Jets.root',
                               'data_pp/MC8TeV_W3Jets.root',
                               ],
                               '#4575b4',
                               [11,13],#[11,13],
                               False )
    if data=='Data8TeV_pp':
        dataset=OrderedDict()
        dataset['l+jets'] = ( ['data_pp/Data8TeV_SingleMu2012A.root'
                                #'data_pp/Data8TeV_SingleMu2012B.root',
                                #'data_pp/Data8TeV_SingleMu2012C.root',
                                #'data_pp/Data8TeV_SingleMu2012D.root'
                                ],
                                 '#4575b4',
                                 [11,13],
                                 False )
        dataset['non-iso l+jets'] = ( dataset['l+jets'][0],'#fc8d59',[11000,13000], False )

    if data=='MC8TeV_TTJets':
        dataset=OrderedDict()
        dataset['l+jets'] = ( ['data_pp/MC8TeV_TTJets_MSDecays_172v5.root'],
                                 '#4575b4',
                                 [13], #[11,13]
                                 True )
        #dataset['non-iso l+jets'] = ( dataset['l+jets'][0], '#91bfdb', [11000,13000],False )

        #dataset['non-iso l+jets'] = ( dataset['l+jets'][0], '#91bfdb', [11000,13000],False )

    return dataset
