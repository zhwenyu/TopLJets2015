import os
import sys
import optparse
import ROOT
import commands
import getpass
import pickle
import numpy
import copy
from subprocess import Popen, PIPE, STDOUT

from TopLJets2015.TopAnalysis.Plot import *
from TopLJets2015.TopAnalysis.dataCardTools import *

def createHistoOrAdd(base,new) :
    if base is None :
        base={}
    for thist in new :
        if thist in base.keys():
            base[thist].Add(new[thist])
        else :
            base[thist]=copy.deepcopy(new[thist])
    return base

def allAreIn(a,b) :
    return len(a) == len([x for x in a if x in b])

def noPartIn(a,b) :
    return 0 == len([c for c in b if (a in c or c in a)])

def buildRateUncertainty(varDn,varUp):
    """ returns a uniformized string for the datacard for the treatment of up/down rate uncertainties """

    #if this is below 0.1% neglect it
    if abs(varUp-1)<0.001 and abs(varDn-1)<0.001: return None

    #distinguish same sided from double sided variations
    toReturn=None
    if varUp>1. and varDn>1.:
        satUnc=min(max(varUp,varDn),2.0)
        toReturn='%3.3f'%satUnc
    elif varUp<1. and varDn<1.:
        satUnc=max(min(varUp,varDn),1/2.0)
        toReturn='%3.3f'%satUnc
    else:
        toReturn='%3.3f/%3.3f'%(varDn,varUp)

    #all done here
    return toReturn

def getDistsFromDirIn(url,indir,applyFilter=''):
    """customize getting the distributions for hypothesis testing"""
    fIn=ROOT.TFile.Open(url)
    obs,exp=getDistsFrom(fIn.Get(indir),applyFilter)
    fIn.Close()
    return obs,exp

def getRowFromTH2(tempHist2D,columnName) :
    """projects a given row in a 2D histogram relying on the y-axis label"""
    tempBinNum = tempHist2D.GetYaxis().FindBin(columnName);
    new1DProj  = tempHist2D.ProjectionX(tempHist2D.GetName()+"_"+columnName,
            tempBinNum,
            tempBinNum) #.Clone(tempHist2D.GetName())
    new1DProj.SetDirectory(0)
    return new1DProj


def getDistsForHypoTest(cat,rawSignalList,opt,outDir="",systName="",systIsGen=False):
    """readout distributions from ROOT file and prepare them to be used in the datacard"""

    # do we want systematics?
    systExt = ""
    if systIsGen        : systExt = "_gen"
    elif systName != "" : systExt = "_exp"

    # get main dists for data and backgrounds
    obs,bkg=getDistsFromDirIn(opt.input,'%s_%s_w100_m1725%s'%(cat,opt.dist,systExt))
    smHypo={}
    for rawSignal in rawSignalList:
        smHypo[rawSignal]=bkg[rawSignal]
        del bkg[rawSignal]

    #get main hypo
    _,mainHypo=getDistsFromDirIn(opt.input,'%s_%s_w%.0f_m%s%s'%(cat,opt.dist,opt.mainHypo,opt.tmass,systExt))
    for proc in [proc for proc in mainHypo if proc in bkg]: del mainHypo[proc]

    #get alternative hypo
    _,altHypo=getDistsFromDirIn(opt.input,'%s_%s_w%.0f_m%s%s'%(cat,opt.dist,opt.altHypo,opt.alttmass,systExt))

    # for 2D measurement: load mass template from correct file, with appropriate names
    if len(opt.altHypoFromSim)>0 :
        _,altHypo=getDistsFromDirIn(opt.systInput,'%s_%s_w%.0f_m%s%s'%(cat,opt.dist,opt.altHypo,opt.alttmass,systExt))
        altHypo = {k.replace(opt.altHypoFromSim,""): v
                        for k, v in altHypo.items()
                        if opt.altHypoFromSim in v.GetName()}

    for proc in [proc for proc in altHypo if proc in bkg]: del altHypo[proc]

    #force the yields to be preserved for alternative hypothesis wrt to the SM expectations
    #if systematics were required we convert the TH2 to the corresponding row
    exp={}
    for proc in bkg:
        if systName != "" :
            bkg[proc] = getRowFromTH2(bkg[proc],systName)
        exp[proc]=bkg[proc]
    for proc in smHypo:
        if systName != "" :
            smHypo[proc] = getRowFromTH2(smHypo[proc],systName)
    for proc in mainHypo:
        if systName != "" :
            mainHypo[proc] = getRowFromTH2(mainHypo[proc],systName)
        nbins=smHypo[proc].GetNbinsX()
        sf=smHypo[proc].Integral(0,nbins+1)/mainHypo[proc].Integral(0,nbins+1)
        mainHypo[proc].Scale(sf)
        newProc='%sw%d'%(proc,opt.mainHypo)
        exp[newProc]=mainHypo[proc].Clone(newProc)
        exp[newProc].SetDirectory(0)
    for proc in altHypo:
        if systName != "" :
            altHypo[proc] = getRowFromTH2(altHypo[proc],systName)
        nbins=altHypo[proc].GetNbinsX()
        sf=smHypo[proc].Integral(0,nbins+1)/altHypo[proc].Integral(0,nbins+1)
        altHypo[proc].Scale(sf)
        newProc='%sw%d'%(proc,opt.altHypo)
        if opt.altHypo==opt.mainHypo : newProc+='a'
        exp[newProc]=altHypo[proc].Clone(newProc)
        exp[newProc].SetDirectory(0)



    return obs,exp,bkg,smHypo


"""
prepare the steering script for combine
"""
def doCombineScript(opt,args,outDir,dataCardList):

    altHypoTag=('w%.0f'%opt.altHypo).replace('.','p')
    if opt.altHypo==opt.mainHypo : altHypoTag+='a'

    scriptname='%s/steerHypoTest.sh'%outDir
    script=open(scriptname,'w')
    print 'Starting script',scriptname
    script.write('#\n')
    script.write('# Generated by %s with git hash %s for standard (alternative) hypothesis %.0f (%.0f)\n' % (getpass.getuser(),
                                                                                                               commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1],
                                                                                                               opt.mainHypo,
                                                                                                               opt.altHypo) )
    script.write('### environment setup\n')
    script.write('COMBINE=%s\n'%opt.combine)
    script.write('SCRIPTDIR=`dirname ${0}`\n')
    script.write('cd ${COMBINE}\n')
    script.write('eval `scramv1 r -sh`\n')
    script.write('cd ${SCRIPTDIR}\n')
    script.write('\n')

    script.write('### combine datacard and start workspace\n')
    script.write('combineCards.py %s > datacard.dat\n'%dataCardList)
    script.write('\n')
    script.write('echo "* autoMCStats 4 0 2" >> datacard.dat')
    script.write('\n')

    script.write('### convert to workspace\n')
    script.write('text2workspace.py datacard.dat -P HiggsAnalysis.CombinedLimit.TopHypoTest:twoHypothesisTest -m 172.5 --PO verbose --PO altSignal=%s --PO muFloating -o workspace.root\n'%altHypoTag)
    script.write('\n')

    def writeScanToScript(testStat,script):
        extraName='_'+testStat
        commonOpts="-m 172.5 -M HybridNew --testStat=%s --onlyTestStat --saveToys --saveHybridResult"%(testStat)
        #commonOpts+=" --cminDefaultMinimizerType GSLMultiMin --cminDefaultMinimizerAlgo BFGS2"
        if hasattr(opt,"frzString") and opt.frzString != "" :
            commonOpts += " --freezeParameters %s"%opt.frzString
        if hasattr(opt,"externStr") and opt.externStr != "" :
            commonOpts += " --setParameters %s"%opt.externStr
        script.write("combine %s --singlePoint 0  workspace.root -n scan0n\n"%commonOpts)
        script.write("mv higgsCombinescan0n.HybridNew.mH172.5.123456.root testStat_scan0n%s.root\n"%extraName)
        script.write("combine %s --singlePoint 1  workspace.root -n scan1n\n"%commonOpts)
        script.write("mv higgsCombinescan1n.HybridNew.mH172.5.123456.root testStat_scan1n%s.root\n"%extraName)

    script.write('### SCAN \n')
    script.write('\n')
    for testStat in ['PL']: writeScanToScript(testStat=testStat,script=script)

    commonOpts="-m 172.5 -M FitDiagnostics --saveWorkspace --saveToys --toysFrequentist --minos none --noErrors"
    #commonOpts+=" --cminDefaultMinimizerType GSLMultiMin --cminDefaultMinimizerAlgo BFGS2"
    setFrzStr="x"
    if hasattr(opt,"frzString") and opt.frzString != "" :
        setFrzStr += ","+opt.frzString
    setParStr="x=1,r=1"
    if hasattr(opt,"externStr") and opt.externStr != "" :
        setParStr += ","+opt.externStr
    script.write("combine %s -t %i --redefineSignalPOIs r --rMin 0.9 --rMax 1.1 --setParameters %s --expectSignal=1 --freezeParameters %s workspace.root -n FitToys\n"%(commonOpts,opt.nToys,setParStr,setFrzStr))

    #repeat fits per categories if validation is required
    script.write('\n')
    script.close()

    return scriptname

"""
instantiates one datacard per category
"""
def doDataCards(opt,args):

    # what are our signal processes?
    rawSignalList=opt.signal.split(',')
    ttScenarioList=['tbart']
    mainSignalList,altSignalList=[],[]
    if 'tbart' in rawSignalList:
        ttScenarioList = ['tbartw%d'%h for h in [opt.mainHypo,opt.altHypo]]
        if opt.mainHypo==opt.altHypo: ttScenarioList[1]+='a'
        mainSignalList += [ttScenarioList[0]]
        altSignalList  += [ttScenarioList[1]]
    tWScenarioList=['Singletop']
    if 'Singletop' in rawSignalList:
        tWScenarioList = ['Singletopw%d'%h for h in [opt.mainHypo,opt.altHypo]]
        if opt.mainHypo==opt.altHypo: tWScenarioList[1]+='a'
        mainSignalList += [tWScenarioList[0]]
        altSignalList  += [tWScenarioList[1]]

    # what are our categories?
    mergedCats=[x.split(';') for x in opt.mergeCatsBy.split(',')]

    # what nuisances do we remove/freeze?
    if opt.rmvNuisances != "" :
        rmvNuisances = True
        nuisanceRMV = opt.rmvNuisances.split(',')
    else :
        rmvNuisances = False

    if opt.frzNuisances != "" :
        frzNuisances = True
        nuisanceFRZ = opt.frzNuisances.split(',')
    else :
        frzNuisances = False


    #define RATE systematics : syst,val,pdf,whiteList,blackList  (val can be a list of values [-var,+var])
    rateSysts=[
          ('lumi_13TeV',       1.025,    'lnN',    [],                  ['DY','W']),
          ('DYnorm',           1.30,     'lnN',    ['DY'],              []),
          ('Wnorm_th',         1.50,     'lnN',    ['W'],               []),
          ('tWnorm_th',        1.15,     'lnN',    tWScenarioList,      []),
          ('VVnorm_th',        1.20,     'lnN',    ['Multiboson'],      []),
          ('tbartVnorm_th',    1.30,     'lnN',    ['tbartV'],          []),
    ]

    #define the SHAPE systematics from weighting, varying object scales, efficiencies, etc.
    # syst,weightList,whiteList,blackList,shapeTreatement=0 (none), 1 (shape only), 2 (factorizeRate),nsigma
    # a - in front of the process name in the black list will exclude rate uncertainties
    weightingSysts=[
        ('ees',            ['ees'],                                    [],             ['DY','-W'], 2, 1.0),
        ('mes',            ['mes'],                                    [],             ['DY','-W'], 2, 1.0),
        ('jer',            ['jer'],                                    [],             ['DY','-W'], 2, 1.0),
        ('trig',           ['trig'],                                   [],             ['DY','-W'], 2, 1.0),
        ('sel_E',          ['esel'],                                   [],             ['DY','-W'], 2, 1.0),
        ('sel_M',          ['msel'],                                   [],             ['DY','-W'], 2, 1.0),
        ('ltag',           ['ltag'],                                   [],             ['DY','-W'], 2, 1.0),
        ('btag',           ['btag'],                                   [],             ['DY','-W'], 2, 1.0),
        ('bfrag',          ['bfrag'],                                  [],             ['DY','-W'], 2, 1.0),
        ('semilep',        ['semilep'],                                [],             ['DY','-W'], 2, 1.0),
        ('pu',             ['pu'],                                     [],             ['DY','-W'], 1, 1.0),
        ('tttoppt',        ['toppt'],                                  ttScenarioList, [],          2, 1.0),
        ('ttMEqcdscale',   ['gen%d'%ig for ig in[3,5,6,4,8,10] ],      ttScenarioList, [],          1, 1.0),
        ('ttPDF',          ['gen%d'%(11+ig) for ig in xrange(0,100) ], ttScenarioList, [],          0, 1.0)
        ]
    for ig in xrange(0,29) :
        weightingSysts += [('jes%s'%ig,            ['jes%d'%ig],       [],             ['DY'], 2, 1.0)]

    #define the SHAPE systematics from dedicated samples : syst,{procs,samples}, shapeTreatment (see above) nsigma
    fileShapeSysts = [
        ('UE',             {'tbart':['t#bar{t} UEdn',     't#bar{t} UEup']}          , 2, 1.0 ),
        ('CR',             {'tbart':['t#bar{t} QCDbased', 't#bar{t} gluon move']}    , 2, 1.0 ),
        ('hdamp',          {'tbart':['t#bar{t} hdamp dn', 't#bar{t} hdamp up']}      , 2, 1.0 ),
        ('ISR_tt',         {'tbart':['t#bar{t} isr dn',   't#bar{t} isr up']}        , 2, 1.0 ),
        ('FSR_tt',         {'tbart':['t#bar{t} fsr dn',   't#bar{t} fsr up']}        , 2, 1.0 ),
        ('ISR_st',         {'Singletop':['Single top isr dn', 'Single top isr up']}  , 2, 1.0 ),
        ('FSR_st',         {'Singletop':['Single top fsr dn', 'Single top fsr up']}  , 2, 1.0 ),
        ('tWttInterf',     {'Singletop':   ['Single top DS']}                        , 2, 1.0 ),
        ('tWMEScale',      {'Singletop':   ['Single top me dn', 'Single top me up']} , 2, 1.0 ),
        ]

    if rmvNuisances and "all" not in nuisanceRMV:
        rateSysts      = [a for a in rateSysts      if noPartIn(a[0],nuisanceRMV)]
        weightingSysts = [a for a in weightingSysts if noPartIn(a[0],nuisanceRMV)]
        fileShapeSysts = [a for a in fileShapeSysts if noPartIn(a[0],nuisanceRMV)]
        print "\n"
        print rateSysts
        print "\n"
        print weightingSysts
        print "\n"
        print fileShapeSysts
    elif rmvNuisances and "all" in nuisanceRMV :
        rateSysts=[]
        weightingSysts=[]
        fileShapeSysts=[]

    # really convoluted, but this was the best way, I promise
    if frzNuisances and "all" not in nuisanceFRZ:

        # collect all correct systematic names, including %sRate
        frzString=",".join([",".join([a[0] for a in rateSysts      if not noPartIn(a[0],nuisanceFRZ)]),
                            ",".join([a[0] for a in weightingSysts if not noPartIn(a[0],nuisanceFRZ)]),
                            ",".join([a[0] for a in fileShapeSysts if not noPartIn(a[0],nuisanceFRZ)])])
        frzString+=","
        frzString+=",".join([",".join([a[0]+"Rate" for a in weightingSysts if not noPartIn(a[0],nuisanceFRZ) and a[4]==2]),
                             ",".join([a[0]+"Rate" for a in fileShapeSysts if not noPartIn(a[0],nuisanceFRZ) and a[2]==2])])
        frzString=frzString.replace(',,,',',')
        frzString=frzString.replace(',,',',')
        frzString=frzString[1:] if frzString[0] == "," else frzString
        frzString=frzString[:-1] if frzString[-1] == "," else frzString

        opt.frzString=frzString
        print "\n"
        print frzString
        print "\n"

        print "\n"
        print len(frzString.split(','))

    elif frzNuisances and "all" in nuisanceFRZ :
        frzString = "all"

        opt.frzString=frzString
        print "\n"
        print frzString
        print "\n"


    # prepare output directory
    outDir='%s/hypotest_%.0fvs%.0f_m%svs%s'%(opt.output, opt.mainHypo,opt.altHypo,opt.tmass,opt.alttmass)
    if opt.pseudoData==-1 : outDir += '_data'
    else:
        outDir += '_%.0f'%opt.pseudoData
        outDir += 'pseudodata'
    os.system('mkdir -p %s'%outDir)
    os.system('rm -rf %s/*'%outDir)

    # prepare output ROOT file
    outFile='%s/shapes.root'%outDir
    fOut=ROOT.TFile.Open(outFile,'RECREATE')
    fOut.Close()

    dataCardList=''
    for mcat in mergedCats :
        print "MERGING",','.join(mcat)
        mname='c'+''.join(mcat)

        tobs,texp,tbkg,tsmHypo=None,None,None,None

        # parse the categories to consider
        for cat in opt.cat.split(','):
            if not allAreIn(mcat,cat) : continue
            print "INCLUDING",cat

            #data and nominal shapes
            obs,exp,bkg,smHypo=getDistsForHypoTest(cat,rawSignalList,opt,outDir)

            texp=createHistoOrAdd(texp,exp)
            tbkg=createHistoOrAdd(tbkg,bkg)
            tsmHypo=createHistoOrAdd(tsmHypo,smHypo)

            #recreate data if requested
            if opt.pseudoData!=-1:
                pseudoSignal=None
                print '\t pseudo-data is being generated',
                if len(opt.pseudoDataFromSim) and opt.systInput:
                    print 'injecting signal from',opt.pseudoDataFromSim
                    pseudoDataFromSim=opt.pseudoDataFromSim.replace('_',' ')
                    _,pseudoSignalRaw=getDistsFromDirIn(opt.systInput,'%s_%s_w%.0f_m%s'%(cat,opt.dist,opt.mainHypo,opt.tmass))
                    pseudoSignal={}
                    pseudoSignal['tbart']=[pseudoSignalRaw[x] for x in pseudoSignalRaw if pseudoDataFromSim in x][0]
                elif len(opt.pseudoDataFromWgt):
                    print 'injecting signal from',opt.pseudoDataFromWgt
                    _,pseudoSignal=getDistsFromDirIn(opt.input,'%s%s_%s_w%.0f_m%s'%(opt.pseudoDataFromWgt,cat,opt.dist,opt.mainHypo,opt.tmass),'t#bar{t}')
                    print pseudoSignal,'%s%s_%s_w%.0f_m%s'%(opt.pseudoDataFromWgt,cat,opt.dist,opt.mainHypo,opt.tmass)
                else:
                    print 'injecting signal from weighted',opt.pseudoData
                    _,pseudoSignal=getDistsFromDirIn(opt.input,'%s_%s_w%.0f_m%s'%(cat,opt.dist,opt.pseudoData,opt.tmass))

                print 'Recreating data from'
                obs.Reset('ICE')
                for proc in bkg:
                    print '\t %s %3.0f'%(proc,bkg[proc].Integral())
                    obs.Add(bkg[proc])

                if not (len(opt.pseudoDataFromSim) and opt.systInput) :
                    pseudoSignal={}
                    _,pseudoSignal=getDistsFromDirIn(opt.input,'%s_%s_w%d_m%s'%(cat,opt.dist,opt.pseudoData,opt.tmass))
                for proc in pseudoSignal:
                    if not proc in rawSignalList : continue
                    print "\t\t Including:", proc, pseudoSignal[proc].GetName()
                    nbins=smHypo[proc].GetNbinsX()
                    sf=smHypo[proc].Integral(0,nbins+1)/pseudoSignal[proc].Integral(0,nbins+1)
                    pseudoSignal[proc].Scale(sf)
                    print '\t %s %3.0f (sf=%3.2f)'%(proc,pseudoSignal[proc].Integral(),sf)

                    obs.Add( pseudoSignal[proc] )

                #round up to integers
                for xbin in xrange(0,obs.GetNbinsX()+2): obs.SetBinContent(xbin,int(obs.GetBinContent(xbin)))
                print '\t Total events in pseudo-data %d'%obs.Integral()

                if tobs is None :
                    tobs=obs.Clone()
                else :
                    tobs.Add(obs)
        print texp

        #start the datacard header
        datacardname='%s/datacard_%s.dat'%(outDir,''.join(mcat))
        dataCardList+='%s=%s '%(mname,os.path.basename(datacardname))

        datacard=open(datacardname,'w')
        print 'Starting datacard',datacardname
        datacard.write('#\n')
        datacard.write('# Generated by %s with git hash %s for analysis with combined categories \n' % (getpass.getuser(),
                                                                                          commands.getstatusoutput('git log --pretty=format:\'%h\' -n 1')[1]) )
        datacard.write('#\n')
        datacard.write('imax *\n')
        datacard.write('jmax *\n')
        datacard.write('kmax *\n')
        datacard.write('-'*50+'\n')
        datacard.write('shapes *        * shapes.root %s_%s/$PROCESS %s_%s_$SYSTEMATIC/$PROCESS\n'%(opt.dist,mname,opt.dist,mname))

        #observation
        datacard.write('-'*50+'\n')
        datacard.write('bin a\n')
        datacard.write('observation %3.1f\n' % tobs.Integral())

        #nominal expectations
        datacard.write('-'*50+'\n')
        datacard.write('\t\t\t %16s'%'bin')
        for i in xrange(0,len(texp)): datacard.write('%15s'%'a')
        datacard.write('\n')
        datacard.write('\t\t\t %16s'%'process')

        for sig in mainSignalList: datacard.write('%15s'%sig)
        for sig in altSignalList:  datacard.write('%15s'%sig)
        for proc in texp:
            if proc in mainSignalList+altSignalList : continue
            datacard.write('%15s'%proc)
        datacard.write('\n')
        datacard.write('\t\t\t %16s'%'process')
        procCtr=-len(mainSignalList)-len(altSignalList)+1
        for sig in mainSignalList:
            datacard.write('%15s'%str(procCtr))
            procCtr+=1
        for sig in altSignalList:
            datacard.write('%15s'%str(procCtr))
            procCtr+=1
        for proc in texp:
            if proc in mainSignalList+altSignalList : continue
            datacard.write('%15s'%str(procCtr))
            procCtr+=1
        datacard.write('\n')
        datacard.write('\t\t\t %16s'%'rate')
        for sig in mainSignalList: datacard.write('%15s'%('%3.2f'%(texp[sig].Integral())))
        for sig in altSignalList:
            #if 'Singletop' in sig : sig = 'Singletopw100'
            datacard.write('%15s'%('%3.2f'%(texp[sig].Integral())))
        for proc in texp:
            if proc in mainSignalList+altSignalList : continue
            datacard.write('%15s'%('%3.2f'%(texp[proc].Integral())))
        datacard.write('\n')
        datacard.write('-'*50+'\n')

        #save to nominal to shapes file
        nomShapes=texp.copy()
        nomShapes['data_obs']=tobs
        saveToShapesFile(outFile,nomShapes,('%s_%s'%(opt.dist,mname)),opt.rebin)

        #rate systematics: these are fixed values common to all processes
        print '\t rate systematics',len(rateSysts)
        for syst,val,pdf,whiteList,blackList in rateSysts:
            datacard.write('%32s %8s'%(syst,pdf))
            entryTxt=''
            try:
                entryTxt='%15s'%('%3.3f/%3.3f'%(ROOT.TMath.Max(val[0],0.01),val[1]))
            except:
                entryTxt='%15s'%('%3.3f'%val)
            for sig in mainSignalList:
                if (len(whiteList)==0 and not sig in blackList) or sig in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for sig in altSignalList:
                if (len(whiteList)==0 and not sig in blackList) or sig in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for proc in texp:
                if proc in mainSignalList+altSignalList : continue
                if (len(whiteList)==0 and not proc in blackList) or proc in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')


        #weighting systema     fortics
        print '\t weighting systematics',len(weightingSysts)
        for syst,weightList,whiteList,blackList,shapeTreatment,nsigma in weightingSysts:
            print '\t\t',syst

            tiexpUp,tiexpDn=None,None
            taltExp,taltExpUp,taltExpDn=None,None,None
            # parse the categories to consider
            for cat in opt.cat.split(','):
                if not allAreIn(mcat,cat) : continue

                isGen = any("gen" in twght for twght in weightList)
                isGen = isGen or 'toppt' in weightList

                # jes has annoying formatting different

                #get shapes and adapt them
                iexpUp,iexpDn=None,None
                altExp,altExpUp,altExpDn=None,None,None
                if len(weightList)==1:
                    if 'jes' in weightList[0] :
                        jesNum=weightList[0].replace('jes','')
                        _,iexpUp,_,_=getDistsForHypoTest(cat,rawSignalList,opt,"","jesup_"+jesNum,isGen)
                        _,iexpDn,_,_=getDistsForHypoTest(cat,rawSignalList,opt,"","jesdn_"+jesNum,isGen)
                    else :
                        _,iexpUp,_,_=getDistsForHypoTest(cat,rawSignalList,opt,"",weightList[0]+"up",isGen)
                        _,iexpDn,_,_=getDistsForHypoTest(cat,rawSignalList,opt,"",weightList[0]+"dn",isGen)

                        if syst=='tttoppt':
                            complCat= cat.replace('lowpt','highpt') if 'lowpt' in cat else cat.replace('highpt','lowpt')
                            _,altExp,_,_   = getDistsForHypoTest(complCat,rawSignalList,opt)
                            _,altExpUp,_,_ = getDistsForHypoTest(complCat,rawSignalList,opt,"",weightList[0]+"up",isGen)

                            #reset the down variation to be the nominal one
                            iexpDn=copy.deepcopy(exp)
                            altExpDn=copy.deepcopy(altExp)

                else:

                    #put all the shapes in a 2D histogram
                    iexp2D={}
                    for iw in xrange(0,len(weightList)):
                        w=weightList[iw]
                        _,kexp,_,_=getDistsForHypoTest(cat,rawSignalList,opt,"",w,isGen)
                        for proc in kexp:
                            nbins=kexp[proc].GetNbinsX()
                            if not proc in iexp2D:
                                name =kexp[proc].GetName()+'2D'
                                title=kexp[proc].GetTitle()
                                xmin =kexp[proc].GetXaxis().GetXmin()
                                xmax =kexp[proc].GetXaxis().GetXmax()
                                nReplicas=len(weightList)
                                iexp2D[proc]=ROOT.TH2D(name,title,nbins,xmin,xmax,nReplicas,0,nReplicas)
                                iexp2D[proc].SetDirectory(0)
                            for xbin in xrange(0,nbins+2):
                                iexp2D[proc].SetBinContent(xbin,iw+1,kexp[proc].GetBinContent(xbin))

                    #create the up/down variations
                    iexpUp,iexpDn={},{}
                    for proc in iexp2D:

                        #create the base shape
                        if not proc in iexpUp:
                            tmp=iexp2D[proc].ProjectionX("tmp",1,1)
                            tmp.Reset('ICE')
                            nbinsx=tmp.GetXaxis().GetNbins()
                            xmin=tmp.GetXaxis().GetXmin()
                            xmax=tmp.GetXaxis().GetXmax()
                            iexpUp[proc]=ROOT.TH1F(iexp2D[proc].GetName().replace('2D','up'),proc,nbinsx,xmin,xmax)
                            iexpUp[proc].SetDirectory(0)
                            iexpDn[proc]=ROOT.TH1F(iexp2D[proc].GetName().replace('2D','dn'),proc,nbinsx,xmin,xmax)
                            iexpDn[proc].SetDirectory(0)
                            tmp.Delete()

                        #project each bin shape for the different variations
                        for xbin in xrange(0,iexp2D[proc].GetNbinsX()+2):
                            tmp=iexp2D[proc].ProjectionY("tmp",xbin,xbin)
                            tvals=numpy.zeros(tmp.GetNbinsX())
                            for txbin in xrange(1,tmp.GetNbinsX()+1) : tvals[txbin-1]=tmp.GetBinContent(txbin)

                            #mean and RMS based
                            if 'PDF' in syst:
                                mean=numpy.mean(tvals)
                                rms=numpy.std(tvals)
                                iexpUp[proc].SetBinContent(xbin,mean+rms)
                                iexpDn[proc].SetBinContent(xbin,ROOT.TMath.Max(mean-rms,1.0e-4))

                            #envelope based
                            else:
                                imax=numpy.max(tvals)
                                if iexpUp[proc].GetBinContent(xbin)>0 : imax=ROOT.TMath.Max(iexpUp[proc].GetBinContent(xbin),imax)
                                iexpUp[proc].SetBinContent(xbin,imax)

                                imin=numpy.min(tvals)
                                if iexpDn[proc].GetBinContent(xbin)>0 : imin=ROOT.TMath.Min(iexpDn[proc].GetBinContent(xbin),imin)
                                iexpDn[proc].SetBinContent(xbin,imin)

                            tmp.Delete()


                        #all done, can remove the 2D histo from memory
                        iexp2D[proc].Delete()

                tiexpUp=createHistoOrAdd(tiexpUp,iexpUp)
                tiexpDn=createHistoOrAdd(tiexpDn,iexpDn)
                taltExp=taltExp if altExp is None else createHistoOrAdd(taltExp,altExp)
                taltExpUp=taltExpUp if altExp is None else createHistoOrAdd(taltExpUp,altExpUp)
                taltExpDn=taltExpDn if altExp is None else createHistoOrAdd(taltExpDn,altExpDn)

                for tbd in iexpUp :
                    iexpUp[tbd].Delete()
                for tbd in iexpDn :
                    iexpDn[tbd].Delete()
                if altExp is not None :
                    for tbd in altExp :
                        if altExp[tbd] is not None : altExp[tbd].Delete()
                    for tbd in altExpUp :
                        if altExpUp[tbd] is not None : altExpUp[tbd].Delete()
                    for tbd in altExpDn :
                        if altExpDn[tbd] is not None : altExpDn[tbd].Delete()

            #write the shapes to the ROOT file
            saveToShapesFile(outFile,tiexpUp,('%s_%s_%sUp'%(opt.dist,mname,syst)),opt.rebin)
            saveToShapesFile(outFile,tiexpDn,('%s_%s_%sDown'%(opt.dist,mname,syst)),opt.rebin)

            #fill in the datacard
            datacard.write('%32s %8s'%(syst,'shape' + (opt.shape)))
            entryTxt='%15s'%('%3.3f'%nsigma)
            for sig in mainSignalList:
                if (len(whiteList)==0 and not sig in blackList) or sig in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for sig in altSignalList:
                if (len(whiteList)==0 and not sig in blackList) or sig in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for proc in texp:
                if proc in mainSignalList+altSignalList : continue
                if (len(whiteList)==0 and not proc in blackList) or proc in whiteList:
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

            #create the rate systematics as well
            iRateVars={}
            if shapeTreatment>0:
                for proc in tiexpUp:
                    nbins=tiexpUp[proc].GetNbinsX()

                    n=texp[proc].Integral(0,nbins+2)
                    nUp=tiexpUp[proc].Integral(0,nbins+2)
                    nDn=tiexpDn[proc].Integral(0,nbins+2)

                    #cases where the rate between high/low pt needs to be taken into account
                    upSF,dnSF=1.0,1.0
                    if taltExp and taltExpUp and taltExpDn:
                        ntot   = n+taltExp[proc].Integral(0,nbins+2)
                        ntotUp = nUp+taltExpUp[proc].Integral(0,nbins+2)
                        ntotDn = nDn+taltExpDn[proc].Integral(0,nbins+2)

                        upSF=ntot/ntotUp
                        dnSF=ntot/ntotDn

                    #normalize shapes to nominal expectations
                    if nUp>0: tiexpUp[proc].Scale(upSF*n/nUp)
                    if nDn>0: tiexpDn[proc].Scale(dnSF*n/nDn)

                    #save a rate systematic from the variation of the yields
                    if n==0 : continue
                    nvarUp=max(nUp/n,nDn/n)
                    nvarDn=min(nUp/n,nDn/n)

                    iRateUnc=buildRateUncertainty(nvarDn,nvarUp)
                    if iRateUnc: iRateVars[proc]=iRateUnc

            #write the rate systematics as well
            if shapeTreatment!=2: continue
            if len(iRateVars)==0: continue
            datacard.write('%32s %8s'%(syst+'Rate','lnN'))
            for sig in mainSignalList:
                if sig in iRateVars and ((len(whiteList)==0 and not sig in blackList and not '-'+sig in blackList) or sig in whiteList):
                    datacard.write('%15s'%iRateVars[sig])
                else:
                    datacard.write('%15s'%'-')
            for sig in altSignalList if opt.useAltRateUncs else mainSignalList:
                if sig in iRateVars and ((len(whiteList)==0 and not sig in blackList and not '-'+sig in blackList) or sig in whiteList):
                    datacard.write('%15s'%iRateVars[sig])
                else:
                    datacard.write('%15s'%'-')
            for proc in texp:
                if proc in mainSignalList+altSignalList : continue
                if proc in iRateVars and ((len(whiteList)==0 and not proc in blackList and not '-'+proc in blackList) or proc in whiteList):
                    datacard.write('%15s'%iRateVars[proc])
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

            for tbd in tiexpUp :
                tiexpUp[tbd].Delete()
            for tbd in tiexpDn :
                tiexpDn[tbd].Delete()


        #systematics from dedicated samples
        print '\t simulated systematics',len(fileShapeSysts)
        for syst,procsAndSamples,shapeTreatment,nsigma in fileShapeSysts:
            print '\t\t',syst

            tiexpUp,tiexpDn=None,None
            # parse the categories to consider
            for cat in opt.cat.split(','):
                if not allAreIn(mcat,cat) : continue

                iexpUp,iexpDn={},{}
                for proc in procsAndSamples:
                    samples=procsAndSamples[proc]

                    hyposToGet=[((opt.mainHypo,opt.tmass) if proc in rawSignalList else (100.0,"1725"))]
                    isSignal=False
                    if proc in rawSignalList:
                        isSignal=True
                        hyposToGet.append((opt.altHypo,opt.alttmass))

                    #use "SM" for non signal
                    if not isSignal:
                        hyposToGet=[(100,"1725")]

                    jexpDn,jexpUp=None,None
                    iInHypos=0
                    for hypo,temptmass in hyposToGet:
                        isAltButMain = iInHypos != 0 and opt.mainHypo == opt.altHypo

                        if len(samples)==2:
                            _,jexpDn=getDistsFromDirIn(opt.systInput,'%s_%s_w%d_m%s'%(cat,opt.dist,hypo,temptmass),samples[0])
                            _,jexpUp=getDistsFromDirIn(opt.systInput,'%s_%s_w%d_m%s'%(cat,opt.dist,hypo,temptmass),samples[1])
                        else:
                            _,jexpUp=getDistsFromDirIn(opt.systInput,'%s_%s_w%d_m%s'%(cat,opt.dist,hypo,temptmass),samples[0])

                        newProc=proc
                        if isSignal: newProc=('%sw%d%s'%(proc,hypo,('a' if isAltButMain else '')))
                        jexpUp.values()[0].SetName(newProc)
                        iexpUp[newProc]=jexpUp.values()[0]

                        #if down variation is not found, mirror it
                        try:
                            jexpDn.values()[0].SetName(newProc)
                            iexpDn[newProc]=jexpDn.values()[0]
                        except:
                            idnHisto=jexpUp.values()[0].Clone()
                            idnHisto.SetDirectory(0)
                            for xbin in xrange(0,idnHisto.GetNbinsX()+2):
                                nomVal=exp[newProc].GetBinContent(xbin)
                                newVal=idnHisto.GetBinContent(xbin)
                                diff=ROOT.TMath.Abs(newVal-nomVal)
                                #if 'tWttInterf' not in syst :
                                #    if newVal>nomVal: nomVal-= ROOT.TMath.Max(diff,1e-4)
                                #    else: nomVal+=diff
                                idnHisto.SetBinContent(xbin,nomVal)
                            iexpDn[newProc]=idnHisto

                        iInHypos+=1

                tiexpUp=createHistoOrAdd(tiexpUp,iexpUp)
                tiexpDn=createHistoOrAdd(tiexpDn,iexpDn)

            #write the shapes to the ROOT file
            saveToShapesFile(outFile,tiexpUp,('%s_%s_%sUp'%(opt.dist,mname,syst)),opt.rebin)
            saveToShapesFile(outFile,tiexpDn,('%s_%s_%sDown'%(opt.dist,mname,syst)),opt.rebin)

            #fill in the datacard
            datacard.write('%32s %8s'%(syst,'shape' + opt.shape))

            for sig in mainSignalList:
                if sig in tiexpUp:
                    entryTxt='%15s'%('%3.3f'%(nsigma if not isinstance(nsigma,dict) else nsigma[sig.split('w')[0]]))
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for sig in altSignalList:
                if 'w%.0fa'%opt.mainHypo in sig :
                    sig=sig.replace('w%.0fa'%opt.mainHypo,'w%.0f'%opt.mainHypo)
                if sig in tiexpUp:
                    entryTxt='%15s'%('%3.3f'%(nsigma if not isinstance(nsigma,dict) else nsigma[sig.split('w')[0]]))
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            for proc in texp:
                if proc in mainSignalList+altSignalList : continue
                if proc in tiexpUp:
                    entryTxt='%15s'%('%3.3f'%(nsigma if not isinstance(nsigma,dict) else nsigma[proc.split('w')[0]]))
                    datacard.write(entryTxt)
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

            #check the shapes
            iRateVars={}
            if shapeTreatment>0:
                for proc in tiexpUp:
                    nbins=tiexpUp[proc].GetNbinsX()

                    #normalize shapes to nominal expectations
                    n=texp[proc].Integral(0,nbins+2)
                    nUp=tiexpUp[proc].Integral(0,nbins+2)
                    if nUp>0: tiexpUp[proc].Scale(n/nUp)
                    nDn=tiexpDn[proc].Integral(0,nbins+2)
                    if nDn>0: tiexpDn[proc].Scale(n/nDn)

                    #save a rate systematic from the variation of the yields
                    if n==0 : continue
                    nvarUp=max(nUp/n,nDn/n) # ROOT.TMath.Abs(1-nUp/n)
                    nvarDn=min(nUp/n,nDn/n) # ROOT.TMath.Abs(1-nDn/n)

                    iRateUnc=buildRateUncertainty(nvarDn,nvarUp)
                    if iRateUnc: iRateVars[proc]=iRateUnc

            #write the rate systematics as well
            if shapeTreatment!=2: continue
            if len(iRateVars)==0 : continue
            datacard.write('%32s %8s'%(syst+'Rate','lnN'))
            for sig in mainSignalList:
                if sig in iRateVars :
                    datacard.write('%15s'%iRateVars[sig])
                else:
                    datacard.write('%15s'%'-')
            for sig in altSignalList if opt.useAltRateUncs else mainSignalList:
                if sig in iRateVars :
                    datacard.write('%15s'%iRateVars[sig])
                else:
                    datacard.write('%15s'%'-')
            for proc in texp:
                if proc in mainSignalList+altSignalList : continue
                if proc in iRateVars :
                    datacard.write('%15s'%iRateVars[proc])
                else:
                    datacard.write('%15s'%'-')
            datacard.write('\n')

            for tbd in tiexpUp :
                tiexpUp[tbd].Delete()
            for tbd in tiexpDn :
                tiexpDn[tbd].Delete()

        print '\t ended datacard generation'
        datacard.close()

    return outDir,dataCardList

"""
steer the script
"""
def main():

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option(      '--combine',            dest='combine',            help='CMSSW_BASE for combine installation',         default=None,        type='string')
    parser.add_option('-i', '--input',              dest='input',              help='input plotter',                               default=None,        type='string')
    parser.add_option(      '--systInput',          dest='systInput',          help='input plotter for systs from alt samples',    default=None,        type='string')
    parser.add_option('-d', '--dist',               dest='dist',               help='distribution',                                default='minmlb',    type='string')
    parser.add_option(      '--shape',              dest='shape',              help='append to shape',                             default='',          type='string')
    parser.add_option(      '--nToys',              dest='nToys',              help='toys to through for CLs',                     default=2000,        type=int)
    parser.add_option(      '--rebin',              dest='rebin',       help='histogram rebin factor',                             default=0,             type=int)
    parser.add_option(      '--pseudoData',         dest='pseudoData',         help='pseudo data to use (-1=real data)',           default=100,         type=float)
    parser.add_option(      '--useAltRateUncs',     dest='useAltRateUncs',     help='use rate uncertainties specific to alt. hypothesis', default=False,       action='store_true')
    parser.add_option(      '--replaceDYshape',     dest='replaceDYshape',     help='use DY shape from syst file',                 default=False,       action='store_true')
    parser.add_option(      '--pseudoDataFromSim',  dest='pseudoDataFromSim',  help='pseudo data from dedicated sample',           default='',          type='string')
    parser.add_option(      '--pseudoDataFromWgt',  dest='pseudoDataFromWgt',  help='pseudo data from dedicated sample',           default='',          type='string')
    parser.add_option(      '--mainHypo',           dest='mainHypo',           help='main width hypothesis',                       default=100,         type=float)
    parser.add_option(      '--tmass',              dest='tmass',              help='main mass hypothesis',                        default=400,         type='string')
    parser.add_option(      '--altHypo',            dest='altHypo',            help='alternative width hypothesis',                default="1725",      type=float)
    parser.add_option(      '--alttmass',           dest='alttmass',           help='alternative mass hypothesis',                 default="1725",      type='string')
    parser.add_option(      '--altHypoFromSim',     dest='altHypoFromSim',     help='alternative hypothesis to take from systs',   default="",          type='string')
    parser.add_option('-s', '--signal',             dest='signal',             help='signal (csv)',                                default='tbart,Singletop',  type='string')
    parser.add_option(      '--removeNuisances',    dest='rmvNuisances',       help='nuisance group to remove (csv)',              default='',  type='string')
    parser.add_option(      '--freezeNuisances',    dest='frzNuisances',       help='nuisance group to freeze (csv)',              default='',  type='string')
    parser.add_option(      '--externNuisances',    dest='externStr',          help='setParameters string',                        default='',  type='string')
    parser.add_option('-c', '--cat',                dest='cat',                help='categories (csv)',
                      default='EE1blowpt,EE2blowpt,EE1bhighpt,EE2bhighpt,EM1blowpt,EM2blowpt,EM1bhighpt,EM2bhighpt,MM1blowpt,MM2blowpt,MM1bhighpt,MM2bhighpt',
                      type='string')
    parser.add_option(      '--mergeCatsBy',        dest='mergeCatsBy',        help='full list of merged categories (colonsv,csv)',
                      default='EE1blowpt,EE2blowpt,EE1bhighpt,EE2bhighpt,EM1blowpt,EM2blowpt,EM1bhighpt,EM2bhighpt,MM1blowpt,MM2blowpt,MM1bhighpt,MM2bhighpt',
                      type='string')
    parser.add_option('-o', '--output',             dest='output',             help='output directory',                            default='datacards', type='string')
    parser.add_option(      '--lumix',              dest='lumix',       help='nuisance group to remove (csv)',              default=1.,  type=float)
    (opt, args) = parser.parse_args()

    outDir,dataCardList=doDataCards(opt,args)
    scriptname=doCombineScript(opt,args,outDir,dataCardList)
    print 'Running statistical analysis'
    runCombine=Popen(['sh',scriptname],stdout=PIPE,stderr=STDOUT)
    runCombine.communicate()

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
