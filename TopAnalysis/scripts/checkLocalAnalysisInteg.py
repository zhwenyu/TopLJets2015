import ROOT
import sys
import os

def modification_date(f):
    t = os.path.getmtime(f)
    return datetime.datetime.fromtimestamp(t)

def checkIntegrity(f,tname):

    if f.find('/store/cms')==0: 
        f='root://eoscms//'+f

    errorCode=None
    nentries=None

    #try to open the file and do a basic check
    try:
        inF=ROOT.TFile.Open(f)
        keys=inF.GetListOfKeys()
        if inF.IsZombie():
            errorCode=1
        if keys.GetSize()==0:
            errorCode=2
    except:
        errorCode=0
        return (errorCode,nentries)
        
    #read tree if available
    if tname:
        try:
            t=inF.Get(tname)
            nentries=t.GetEntriesFast()
        except: 
            errorCode=3

    #close file
    inF.Close()
    
    tstamp=None
    try:
        tsamp=modification_date(f.replace('root://eoscms//','/eos/cms'))
    except:
        pass

    return (errorCode,nentries,tstamp)




#read files
f=open(os.path.join(sys.argv[1],'checkIntegList.dat'))
files = [line.rstrip('\n').split() for line in f]
f.close()
print len(files),'files to check'

itname=sys.argv[2] if len(sys.argv)>2 else None
otname=sys.argv[3] if len(sys.argv)>3 else None

with open('localanalysis_integ_report.dat','w') as r:
    for i,o in files:
        ires=checkIntegrity(i,itname)
        ores=checkIntegrity(o,otname)
        msg='IS OK' if ires[0] is None and ores[0] is None else 'NOTOK'
        ijob='/'.join(i.split('/')[-2:])
        r.write('[{0:6s}] job with input {1} input={2} output={3}\n'.format(msg,ijob,ires,ores))

print 'Report in localanalysis_integ_report.dat'

