import ROOT
import sys
import os
import re

def sumUpAllAnglesFor(boson,inDir):

    #list all files for this boson
    flist=[os.path.join(inDir,f) for f in os.listdir(inDir) if boson in f and '.root' in f]

    plots={}
    print 'Merging plots in',flist
    for f in flist:
        inF=ROOT.TFile.Open(f)
        for k in inF.GetListOfKeys():
            newName=k.GetName()
            angle=int(re.search('_a(\d+)',newName).group(1))
            newName=newName.replace('a%d_'%angle,'')
            if newName in plots:
                plots[newName].Add( k.ReadObj() )
            else:
                plots[newName]=k.ReadObj().Clone(newName)
                plots[newName].SetDirectory(0)
        inF.Close()

    #dump to new file
    os.system('mkdir -p %s/inclusive'%inDir)
    fOut=ROOT.TFile.Open(os.path.join(inDir,'inclusive/shapes_%s.root'%boson),'RECREATE')
    for k in plots:
        plots[k].SetDirectory(fOut)
        plots[k].Write()
    fOut.Close()

def generateDataCard(boson,inDir):

    #list all datacards for this boson and one xangle
    bosonstr='zmm'
    if boson==121: bosonstr='zee'
    if boson==22 : bosonstr='g'
    flist=[os.path.join(inDir,f) for f in os.listdir(inDir) if bosonstr in f and '.dat' in f and '_a120' in f]
    print 'Converting to inclusive datacards from',flist

    #sed all instances of the angle and we should be good to go
    combstr=''
    for ictr,dc in list(enumerate(flist)):
        newDC=os.path.basename(dc)
        newDC=inDir+'/inclusive/'+newDC.replace('_a120','')
        os.system("sed 's/_a120//g' {0} > {1}".format(dc,newDC))
        combstr+='cat%d=%s '%(ictr,newDC)
    
    #combine datacard
    combine='combineCards.py {0} > {1}/inclusive/{2}_datacard.dat'.format(combstr,inDir,bosonstr)
    os.system(combine)

for boson in [11*11,13*13,22]:
    sumUpAllAnglesFor(str(boson),sys.argv[1])
    generateDataCard(boson,sys.argv[1])
