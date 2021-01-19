import ROOT
import os
import sys
from itertools import product
from runExclusiveAnalysis import getTracksPerRomanPot,buildDiProton
from TopLJets2015.TopAnalysis.HistoTool import *
from TopLJets2015.TopAnalysis.myProgressBar import *
from TopLJets2015.TopAnalysis.Plot import *

def buildControlPlots(args):
    
    """loop over data events and fill some basic control plots based on Z->mumu pT<10 GeV data"""

    tag,ch=args
    evEra='2017'+tag
    baseDir='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/final/2017_unblind_multi/Chunks/'
    tree=ROOT.TChain('tree')
    for f in os.listdir(baseDir):
        if ch==13*13 and not 'Data13TeV_2017%s_DoubleMuon'%tag in f : continue
        if ch==11*11 and not 'Data13TeV_2017%s_DoubleEG'%tag in f : continue
        tree.AddFile(os.path.join(baseDir,f))

    nentries=tree.GetEntries()
    if nentries==0: return
    print 'Analysing',nentries,'events for tag',tag,'ch=',ch

    #book histograms
    ht=HistoTool()
    ht.add(ROOT.TH1F('n',       ';Proton multiplicity;PDF',                         5,0,5))
    ht.add(ROOT.TH1F('csi',     ';#xi;PDF',                                         50,0,0.3))
    ht.add(ROOT.TH1F('mpp',     ';m_{pp} [GeV];PDF',                                100,0,2500))
    ht.add(ROOT.TH1F('mmass',   ';Missing mass [GeV];Events',                       100,-1000,3000))
    ht.add(ROOT.TH1F('xangle',  ';Beam crossing angle [#murad];Events',             4,120,160))
    ht.add(ROOT.TH1F('nvtx',    ';Vertex multiplicity;Events',                      50,0,50))
    ht.add(ROOT.TH1F('ptll',    ';Transverse momentum [GeV];Events',                20,0,10))
    ht.add(ROOT.TH1F('met',     ';Missing transverse energy [GeV];Events',          20,0,200))
    ht.add(ROOT.TH1F('dphimetz',';#Delta#phi[E_{T}^{miss},p_{T}(ll)] [rad];Events', 20,0,3.15))
    ht.add(ROOT.TH1F('nch',     ';Charged particle multiplicity;Events',            20,0,100))
    ht.add(ROOT.TH1F('rue',     ';p_{T}(vtx)/p_{T}(ll)-1;Events',                   40,-1,1))

    for i in range(nentries):

        tree.GetEntry(i)

        if i%1000==0 : drawProgressBar(float(i)/float(nentries))
        if not tree.isZ : continue
        if abs(tree.l1id*tree.l2id)!=ch: continue
        if tree.bosonpt>10 : continue

        evRun=tree.run
        xangle=tree.beamXangle

        #decode variables of interest
        boson=ROOT.TLorentzVector(0,0,0,0)
        boson.SetPtEtaPhiM(tree.bosonpt,tree.bosoneta,tree.bosonphi,tree.mboson)
        tkPos,tkNeg=getTracksPerRomanPot(tree,evEra,evRun,xangle)
        rue=tree.sumPVChPt/boson.Pt()-1
        dphimetz=abs(ROOT.TVector2.Phi_mpi_pi(tree.met_phi-boson.Phi()))
        xangle=tree.beamXangle

        cats=['inc']        
        for ialgo in range(3):
            algoCat=['multi']
            if ialgo==1: algoCat=['px']
            if ialgo==2: algoCat=['strip']        

            #individual RP plots
            for csiList,side in [(tkPos[ialgo],'pos'),(tkNeg[ialgo],'neg')]:
                ht.fill((len(csiList),1),'n',algoCat,side)
                if len(csiList)==0: continue
                ht.fill((csiList[0],1),'csi',algoCat,side)
                for x in csiList:
                    ht.fill((x,1),'csi',algoCat,side+'_inc')

            #combined PPS variables
            passSel=True if len(tkPos[ialgo])>0 and len(tkNeg[ialgo])>0 else False
            if not passSel: continue
            cats += algoCat

            csi0,csi1=tkPos[ialgo][0],tkNeg[ialgo][0]
            pp=buildDiProton(csi0,csi1)
            ht.fill( (pp.M(),1), 'mpp', algoCat)

            mmass=(pp-boson).M()
            ht.fill( (mmass,1),  'mmass', algoCat)

        #central, global variables
        ht.fill((tree.nvtx,1),        'nvtx',     cats + ['a%d'%xangle])
        ht.fill((xangle,1),           'xangle',   cats)
        ht.fill((boson.Pt(),1),       'ptll',     cats)
        ht.fill((tree.nchPV,1),       'nch',      cats)
        ht.fill((rue,1),              'rue',      cats)
        ht.fill((tree.met_pt,1),      'met',      cats)
        ht.fill((dphimetz,1),         'dphimetz', cats)

    outURL='RPcontrol_era%s_ch%d.root'%(tag,ch)
    ht.writeToFile(outURL)
    print 'Results can be found in',outURL

def drawPlots():

    """draw final plots comparing distributions in data"""

    def fillHistos(fIn,histos,era,ch):

        eraKey=(era,ch)
        for k in fIn.GetListOfKeys():
            kname =k.GetName()
            hname =kname.split('_')[0]
            tag   =''.join(kname.split('_')[1:])
            key=(hname,tag)
            if not key in histos:
                histos[key]={}
            histos[key][eraKey]=k.ReadObj().Clone(kname+'%s_%d'%(era,ch))
            histos[key][eraKey].SetDirectory(0)
            histos[key][eraKey].GetYaxis().SetTitle('PDF')

        return histos

    histos={}
    for era,ch in list(product( list('BCDEF'), [11*11,13*13] )):
        fIn=ROOT.TFile.Open('RPcontrol_era%s_ch%d.root'%(era,ch))
        fillHistos(fIn,histos,era,ch)
        fIn.Close()


    perAlgo=['multi','px','strip']
    perAlgoPlusInc=perAlgo+['']+['a120','a130','a140','a150']
    perSide=product(perAlgo,['pos','posinc','neg','neginc'])
    
    colors=["#fdc086","#7fc97f","#DCDCDC",633,"#386cb0",9]
    for d,tags in [('n',        perSide),
                   ('csi',      perSide),
                   ('mpp',      perAlgo),
                   ('mmass',    perAlgo),
                   ('nvtx',     perAlgoPlusInc),
                   ('xangle',   perAlgoPlusInc),
                   ('ptll',     perAlgoPlusInc),
                   ('nch',      perAlgoPlusInc),
                   ('rue',      perAlgoPlusInc),
                   ('met',      perAlgoPlusInc),
                   ('dphimetz', perAlgoPlusInc),
                   ]:

        hglobal=None
        hperAlgo={}
        hperXangle={}
        for tag in tags:
            hinc=None
            hperEra={}
            hperCh={}

            key=(d,''.join(tag))
            if not key in histos: continue

            for eraKey,h in histos[key].items():

                era,ch=eraKey
                    
                #global inclusive histogram
                if len(tag)==0:
                    if hglobal is None:
                        hglobal=h.Clone(d+'_global')
                        hglobal.SetDirectory(0)
                    else:
                        hglobal.Add(h)

                #inclusive histogram (this tag only)
                if hinc is None:
                    hinc=h.Clone(d+'_inc')
                    hinc.SetDirectory(0)
                else:
                    hinc.Add(h)

                #per angle histos
                if len(tag)==4 and tag[0]=='a':
                    xangle=tag[1:]
                    if not xangle in hperXangle:
                        hperXangle[xangle]=h.Clone(d+'_'+xangle)
                        hperXangle[xangle].SetDirectory(0)
                    else:
                        hperXangle[xangle].Add(h)

                #per era histogram
                if  not era in hperEra:
                    hperEra[era]=h.Clone('%s_%s'%(d,era))
                    hperEra[era].SetDirectory(0)
                else:
                    hperEra[era].Add(h)

                #per ch histogram
                if not ch in hperCh:
                    hperCh[ch]=h.Clone('%s_%d'%(d,ch))
                    hperCh[ch].SetDirectory(0)
                else:
                    hperCh[ch].Add(h)

                #per algo histogram
                algo=tag
                algoExtra=None
                if d=='csi' and not 'inc' in tag[1]:  
                    algo=tag[0]
                    algoExtra =  '_pos' if 'pos' in tag[1] else '_neg'
                if isinstance(algo,str) and (algo in perAlgo+['']):
                    algoKey=algo
                    if algoExtra: algoKey+=algoExtra
                    if not algoKey in hperAlgo:
                        hperAlgo[algoKey]=h.Clone('%s_%s'%(d,algoKey))
                        hperAlgo[algoKey].SetDirectory(0)
                    else:
                        hperAlgo[algoKey].Add(h)


            if hinc is None: continue
            hinc.Scale(1./hinc.Integral())

            #per era
            p=Plot('%s_%s_perEra'%key,com='13 TeV')
            p.savelog=True
            p.range=[1e-3,1]
            p.doPoissonErrorBars=False
            for i,era in enumerate(hperEra.keys()):
                hperEra[era].Scale(1./hperEra[era].Integral())
                p.add(hperEra[era], title=era, color=colors[i], isData=False, spImpose=False, isSyst=False)
            p.add(hinc.Clone(), title='total', color=1,         isData=False, spImpose=True,  isSyst=False)
            p.show(outDir='./', lumi=37500, noStack=True)

            #per ch
            p=Plot('%s_%s_perCh'%key,com='13 TeV')
            p.savelog=True
            p.range=[1e-3,1]
            p.doPoissonErrorBars=False
            for i,ch in enumerate(hperCh.keys()):
                hperCh[ch].Scale(1./hperCh[ch].Integral())
                p.add(hperCh[ch], title='ee' if ch==11*11 else '#mu#mu', color=colors[i], isData=False, spImpose=False, isSyst=False)
            p.add(hinc.Clone(), title='total', color=1,         isData=False, spImpose=True,  isSyst=False)
            p.show(outDir='./', lumi=37500, noStack=True)

        #per xangle
        p=Plot('%s_perxangle'%d,com='13 TeV')
        p.savelog=True
        p.range=[1e-3,1]
        p.doPoissonErrorBars=False
        for i,xangle in enumerate(['120','130','140','150']):
            if not xangle in hperXangle: continue
            hperXangle[xangle].Scale(1./hperXangle[xangle].Integral())
            spImpose=True if xangle=='' else False
            title='%s#murad'%xangle
            ci=colors[i]
            p.add(hperXangle[xangle], title=title, color=ci, isData=False, spImpose=spImpose, isSyst=False)
        if hglobal:
            hglobal.Scale(1./hglobal.Integral())
            p.add(hglobal.Clone(), title='total',         color=1,  isData=False, spImpose=True,     isSyst=False)
        try:
            p.show(outDir='./', lumi=37500, noStack=True)
        except Exception as e:
            print d,e
            
        #per algo
        p=Plot('%s_peralgo'%d,com='13 TeV')
        p.savelog=True
        p.range=[1e-3,1]
        p.doPoissonErrorBars=False
        for i,algo in enumerate(hperAlgo.keys()):
            hperAlgo[algo].Scale(1./hperAlgo[algo].Integral())
            spImpose=True if algo=='' else False
            title='inclusive' if algo=='' else algo
            title=title.replace('_pos',' (+)')
            title=title.replace('_neg',' (-)')
            ci=colors[i] if algo!='' else 1
            p.add(hperAlgo[algo], title=title, color=ci, isData=False, spImpose=spImpose, isSyst=False)
        p.show(outDir='./', lumi=37500, noStack=True)


def main() :

    if sys.argv[1]=='fill':
        import multiprocessing as MP

        pool = MP.Pool(8)
        task_list = list(product( list('BCDEF'), [11*11,13*13] ))
        pool.map(buildControlPlots,task_list)

    else:
        ROOT.gROOT.SetBatch(True)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetPaintTextFormat("4.2f");
        drawPlots()


if __name__ == "__main__":
    sys.exit(main())
