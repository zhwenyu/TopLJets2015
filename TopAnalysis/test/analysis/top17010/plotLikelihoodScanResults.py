import ROOT
import os
import sys
import optparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
import re

def getTheoryPrediction(x=np.arange(169,176,0.1)):
    
    """top width versus top mass at NLO """

    def nlo_gt(mt):
        GF=1.1663787e-5
        mW=80.385
        aS=0.118
        pi=ROOT.TMath.Pi()
        
        gt  = GF*(mt**3)/(8*pi*ROOT.TMath.Sqrt(2.))
        gt *= (1-(mW/mt)**2)**2
        gt *= (1+2*(mW/mt)**2)
        gt *= (1-(2*aS/(3*pi))*(2*(pi**2)/3-5./2.))

        return gt

    return x,np.array( [nlo_gt(m) for m in x] )
    




def getScanPoint(inDir,fitTag):

    """read the fit result and return the likelihood value together with the corresponding mtop,width values"""

    #decode mass and width
    tag=os.path.basename(inDir)
    mt,gt=172.5,1.31
    if 'scenario' in tag:
        flag=int(re.findall(r'\d+',tag)[0])
        mask=int('0xffff',16)
        gt = (flag&mask)*0.01+0.7
        mt = (((flag>>16)&mask))*0.25+169        

    #read nll from fit
    nll=None
    try:
        url=os.path.join(inDir,'fitresults%s.root'%fitTag)
        if os.path.isfile(url):
            inF=ROOT.TFile.Open(url)
            tree=inF.Get('fitresults')
            tree.GetEntry(0)
            nll=tree.nllvalfull
    except Exception as e:
        print e

    return mt,gt,nll
        

def profilePOI(data,outdir,axis=0):

    """ profiles in x and y the POI """

    #raw values
    x=data[:,axis]
    xtit='$m_{t}$ [GeV]' if axis==1 else '$\Gamma_{t}$ [GeV]'
    ytit='$m_{t}$ [GeV]' if axis==0 else '$\Gamma_{t}$ [GeV]'

    ictr=0
    for xi in np.unique(x):
        rdata=data[data[:,axis]==xi]
        y=rdata[:,0 if axis==1 else 1]
        z=rdata[:,2]
        pcoeff=np.polyfit(y,z,4)
        p=np.poly1d(pcoeff)


        bounds = [min(y),max(y)]
        crit_points = [px for px in p.deriv().r if px.imag == 0 and bounds[0] < px.real < bounds[1]]
        print xi,crit_points

        plt.clf()
        fig, ax = plt.subplots()
        yp = np.linspace(bounds[0],bounds[1], 100)
        plt.plot(y, z, '.',label='scan points')
        plt.plot(yp, p(yp), '-', label='interpolation')
        plt.xlabel(xtit)
        #plt.ylim(0.,20.0)
        ax.text(0,1.02,'CMS preliminary', transform=ax.transAxes, fontsize=16)
        ax.text(1.0,1.02,r'%s=%3.2f 34.5 fb$^{-1}$ (13 TeV)'%(ytit,xi), transform=ax.transAxes,horizontalalignment='right',fontsize=14)
        ax.legend(framealpha=0.0, fontsize=14, loc='upper left', numpoints=1)        
        plt.savefig(os.path.join(outdir,'nllprofile_%d_%d.png'%(axis,ictr)))
        ictr+=1



def doContour(data,outdir,              
              method='linear',
              levels=[2.30,4.61,9.21],
              levelLabels=['68.3%','90%','99%'],
              linestyles=['-','-','-']):

    """ interpolates the grid to obtain the likelihood contour 
    2 parameter fit levels (see PDG Statistics Table 38.2) """

    #raw values
    x=data[:,0]
    y=data[:,1]
    z=data[:,2]

    #interpolate and find minimum
    xi = np.linspace(169.5, 175.5,100)
    yi = np.linspace(0.7,4.0,100)
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method=method)
    minz=zi.min()
    bestFitIdx = np.where(zi==minz)
    bestX=xi[ bestFitIdx[1][0] ]
    bestY=yi[ bestFitIdx[0][0] ]
    zi=(zi-minz)*2
    z=(z-minz)*2

    fig, ax = plt.subplots()
    cntr=ax.contour(xi, yi, zi, levels=levels, linewidths=1.0, linestyles=linestyles, colors='k')
    cntr.collections[0].set_label('$-2\Delta ln(\lambda)$')

    #add contour level names
    fmt = {}
    for l, s in zip(levels,levelLabels):
        fmt[l] = s
    ax.clabel(cntr, cntr.levels[:], inline=True, fmt=fmt, fontsize=10)

    #add best-fit point
    plt.plot([bestX], [bestY], '+', mew=4, markersize=12, color='k',label='Best fit')

    #add theory prediction
    theory=getTheoryPrediction()
    plt.plot(theory[0],theory[1],'-',color='lightgray',label='Theory NLO')

    plt.xlabel('$m_{t}$ [GeV]', fontsize=16)
    plt.ylabel('$\Gamma_{t}$ [GeV]', fontsize=16)
    plt.ylim(0.7,4.0)
    plt.xlim(169.5,175.5)
    ax.text(0,1.02,'CMS preliminary', transform=ax.transAxes, fontsize=16)
    ax.text(1.0,1.02,r'34.5 fb$^{-1}$ (13 TeV)', transform=ax.transAxes,horizontalalignment='right',fontsize=14)
    ax.legend(framealpha=0.0, fontsize=14, loc='upper left', numpoints=1)

    for ext in ['png','pdf']:
        plt.savefig(os.path.join(outdir,'nllcontour.%s'%ext))


def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',          
                      dest='input',       
                      help='input directory [%default]',  
                      default='/eos/cms/store/cmst3/group/top/TOP17010/0c522df/fits/em_inc',
                      type='string')
    parser.add_option('-o', '--out',          
                      dest='outdir',
                      help='output directory [%default]',  
                      default='fit_results/em_inc',
                      type='string')
    parser.add_option('-t', '--tag',          
                      dest='fitTag',
                      help='fit tag [%default]',  
                      default='_tbart',
                      type='string')
    (opt, args) = parser.parse_args()


    #build nll scan
    fitres=[]
    for f in os.listdir(opt.input):
        scanRes=getScanPoint(inDir=os.path.join(opt.input,f),fitTag=opt.fitTag)
        if not scanRes[-1]: continue
        fitres.append( scanRes )
    fitres=np.array(fitres)

    #plot the contour interpolating the available points
    os.system('mkdir -p %s'%opt.outdir)
    doContour(fitres,outdir=opt.outdir)
    profilePOI(fitres,outdir=opt.outdir,axis=0)
    profilePOI(fitres,outdir=opt.outdir,axis=1)

if __name__ == "__main__":
    sys.exit(main())
