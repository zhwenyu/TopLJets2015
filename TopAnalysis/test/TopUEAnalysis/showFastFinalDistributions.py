import ROOT
import os
import sys
import optparse
import pickle

from UETools import getPurStab

"""
"""
def readPlotsFrom(args,opt,sliceVar,axisVar):

    outdir=args[0].replace('.root','')
    os.system('mkdir -p %s'%outdir)

    #analysis axes
    analysisaxis=None
    with open(opt.analysisAxis,'r') as cachefile:
        analysisaxis = pickle.load(cachefile)

    c=ROOT.TCanvas('c','c',1000,1000)
    c.Divide(2,2)

    fIn=ROOT.TFile.Open(args[0])
    for obs in ['chmult','sphericity','C','D','aplanarity','chavgpt','chavgpz','chflux','chfluxz']:
        
        hfinal=[]
        for level in [False,True]: 
            obsAxis=analysisaxis[(obs,level)]
            key='%s_%s_%s_None_%s'%(obs,sliceVar,axisVar,level)
            h=fIn.Get(key)
            
            p=c.cd(int(level)+1)
            p.Clear()
            try:
                hfinal.append( h.Clone() )
        
            #hfinal.append( ROOT.TH1F('%s_finalbin'%key,
            #                         '%s_finalbin'%key,
            #                         obsAxis.GetNbins(),
            #                         obsAxis.GetXbins().GetArray()) )
            #
            #for xbin in xrange(1,obsAxis.GetNbins()+1):
            #    val=h.GetBinContent(xbin)
            #    unc=h.GetBinError(xbin)
            #    wid=obsAxis.GetBinWidth(xbin)
            #    hfinal[-1].SetBinContent(xbin,val/wid)
            #    hfinal[-1].SetBinError(xbin,unc/wid)
    
                hfinal[-1].Draw()
            except:
                print 'Failed for',key
                pass

        p=c.cd(3)
        p.Clear()
        key='%s_%s_inc_0_mig'%(obs,sliceVar)
        h2d=fIn.Get(key)
        try:
            hfinal.append( h2d.Clone(key+'_norm') )
            hfinal[-1].Draw('colz')
            for xbin in xrange(1,h2d.GetNbinsX()+1):
                tmp=h2d.ProjectionY('tmp',xbin,xbin)
                total=tmp.Integral(1,tmp.GetNbinsX())
                tmp.Delete()
                if total==0 : continue
                for ybin in xrange(1,h2d.GetNbinsY()+1):
                    val=h2d.GetBinContent(xbin,ybin)
                    unc=h2d.GetBinError(xbin,ybin)
                    hfinal[-1].SetBinContent(xbin,ybin,100.*val/total)
                    hfinal[-1].SetBinError(xbin,ybin,100.*unc/total)
        except:
            pass

        p=c.cd(4)
        p.Clear()
        try:
            h2d.RebinY()
            purGr,stabGr,effGr=getPurStab(h2d)
            stabGr.Draw('al')
            stabGr.GetYaxis().SetRangeUser(0,1)
            purGr.Draw('l')
            purGr.SetLineWidth(2)
            effGr.SetLineColor(ROOT.kGray)
            effGr.Draw('l')
        except:
            pass
        c.cd()
        c.Modified()
        c.Update()
        c.SaveAs('%s/%s_%s_preview.png'%(outdir,obs,sliceVar))
        for h in hfinal: h.Delete()


"""
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--cfg',
                      dest='analysisAxis',
                      help='cfg with axis definitions [%default]',
                      type='string',
                      default='%s/src/TopLJets2015/TopAnalysis/UEanalysis/analysisaxiscfg.pck'%os.environ['CMSSW_BASE'])
    (opt, args) = parser.parse_args()

    for sliceVar in [None] : #,'ptll','ptttbar','nj']:
        for axisVar in ['inc','ptll','ptttbar']:
            readPlotsFrom(args,opt,sliceVar,axisVar)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
