import ROOT
import os
import sys
import numpy as np

systList=['bkgShapeUp',           'bkgShapeDown',
          'sigShapeUp',           'sigShapeDown',
          'bkgShapeSingleDiffUp', 'bkgShapeSingleDiffDown',
          'sigShapeSingleDiffUp', 'sigShapeSingleDiffDown',
          'sigCalibUp',           'sigCalibDown',
          'sigPzModelUp',         'sigPzModelDown']
          

def printUncertainties(flist,m=1000):

    uncVals={}

    for f in flist:
        r=ROOT.TFile.Open(f)
        baseH={}
        for k in r.GetListOfKeys():
            kname=k.GetName()

            if 'data' in kname : continue
            if 'sig_' in kname:
                sigm=int(kname.split('_')[2].replace('m',''))
                if m!=sigm : continue
            
            isSigFid    = True if 'fidsig' in kname and not 'outfidsig' in kname else False
            isSigOutFid = True if 'outfidsig' in kname else False
            sname       = kname.split('_')[-1]
            if not sname in systList: sname='mcstats'
            sname=sname.replace('Up','')
            sname=sname.replace('Down','')

            key=None
            if isSigFid:
                key=('sigfid',sname)
            elif isSigOutFid:
                key=('sigoutfid',sname)
            else:
                key=('bkg',sname)
            if not key in uncVals:
                uncVals[key]=[]

            h=k.ReadObj()
            if key[1]=='mcstats':
                baseH[key[0]]=k.ReadObj()
            else:
                h.Divide(baseH[key[0]])
                
            for xbin in xrange(1,h.GetNbinsX()):
                val=h.GetBinContent(xbin)
                if key[1]=='mcstats':
                    if val<=0.01 :continue
                    uncVals[key].append( h.GetBinError(xbin)/val )
                else:
                    uncVals[key].append( h.GetBinContent(xbin)-1 )

        r.Close()

    for key in uncVals:
        print '%30s %30s \t\t %3.3f~~{\small $]%3.3f,%3.3f[$}'%(key[0],
                                                    key[1],
                                                    np.percentile(uncVals[key], 50),
                                                    np.percentile(uncVals[key], 10),
                                                    np.percentile(uncVals[key], 90))

url=sys.argv[1]

f=[(x.split('_')[1],os.path.join(url,x)) for x in os.listdir(url) if 'shapes_' in x]
chf={}
for ch,r in f:
    if not ch in chf: chf[ch]=[]
    chf[ch].append(r)

for ch in chf:

    print '-'*50
    print ch
    print '-'*50
    printUncertainties(chf[ch])
