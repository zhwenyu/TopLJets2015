import ROOT
import sys
import os
import pickle
import re
from prepareOptimScanCards import OPTIMLIST

def readLimitsFrom(url,getObs=False):

    """parses the r95 limits from the tree"""

    fIn=ROOT.TFile.Open(url)
    try:
        t=fIn.Get('limit')
        vals=[]
        if getObs:
            t.GetEntry(t.GetEntriesFast()-1)
            vals=[t.limit]*5
        else:
            for i in range(5):
                t.GetEntry(i)
                vals.append(t.limit)
        fIn.Close()
    except:
        vals=[999.]*5

    return [vals[2],vals[3]-vals[2],vals[1]-vals[2],vals[4]-vals[2],vals[0]-vals[2]]


def readSignificanceFrom(url):

    """parses the significance value from combine tree and converts it to a p-value"""

    fIn=ROOT.TFile.Open(url)
    try:
        t=fIn.Get('limit')
        t.GetEntry(0)
        vals=[t.limit,ROOT.RooStats.SignificanceToPValue(t.limit)]
        fIn.Close()
    except:
        vals=[0,0.5]

    return vals


def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)

    baseDir=sys.argv[1]

    #scan all results in directory
    results=[]
    toCheck=[]
    for ana in os.listdir(baseDir):

        if 'optim_' in ana: continue
        anaDir=os.path.join(baseDir,ana)
        if not os.path.isdir(anaDir) : continue
        
        limitFiles=[os.path.join(anaDir,x) for x in os.listdir(anaDir) if 'X.obs.AsymptoticLimits' in x]
        for f in limitFiles:
            ch        = re.search('PP([egmz]+)X', f).group(1)
            mass      = int(re.search('.mH(\d+)', f).group(1))
            iresults  = [ana,ch,mass] 
            iresults += readLimitsFrom(f)
            iresults += readSignificanceFrom(f.replace('X.obs.AsymptoticLimits','X.Significance'))
            if iresults[3]<900:
                results.append( iresults )
            else:
                toCheck.append( (ch,mass,anaDir) )

    #save summary in a pandas dataformat
    import pandas as pd
    columns=['ana','channel','mass','r95','drup68','drdn68','drup95','drdn95','sig','pval']
    df=pd.DataFrame(data=results, columns=columns)
    df.to_hdf('%s/summary.h5'%baseDir,key='scan')

    print df.head()
    print 'All results are available in %s/summary.h5'%baseDir

    if len(toCheck)>0:
        print 'Recovering %d missing jobs'%len(toCheck)
        for ch,m,anaDir in toCheck:
            print ch,m,anaDir
            os.system('sh %s/statAnaJob.sh %s %s'%(anaDir,m,ch))


if __name__ == "__main__":
    main()
