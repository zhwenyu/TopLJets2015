import ROOT
import sys
import os

baseDir=sys.argv[1]
optimList=[os.path.join(baseDir,d) for d in os.listdir(baseDir) if 'optim_' in d] 
toCheck=[]
for d in optimList:

    isOK=True
    for ch in [121,169,22]:

        try:
            fIn=ROOT.TFile.Open(os.path.join(d,'shapes_%d.root'%ch))
            size=fIn.GetListOfKeys().GetSize()
            if size==0: isOK=False
            fIn.Close()
        except:
            isOK=False
            pass

        if not isOK: break

    if isOK: continue
    toCheck.append( d.split('_')[-1] )

if len(toCheck)>0:

    resub=os.path.join(baseDir,'zxstatana_scan_resub.sub')
    with open(resub,'w') as f:
        f.write('executable  = %s/optim_$(idx)/optimJob.sh\n'%os.path.abspath(baseDir))
        f.write('output       = zxstatana_scan_resub.out\n')
        f.write('error        = zxstatana_scan_resub.err\n')
        f.write('log          = zxstatana_scan_resub.log\n')
        f.write('+JobFlavour = "tomorrow"\n')
        f.write('request_cpus = 4\n')
        f.write('queue idx in ( %s )\n'%' '.join(toCheck))

    print '%d resubmission jobs can be found in %s'%(len(toCheck),resub)
