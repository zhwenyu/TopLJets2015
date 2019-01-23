import ROOT
import os
import generateWorkspace

table=[]
for l1pt in [30,40,50]:
    for l2pt in [20,30,40,50]:
        for bosonpt in [20,40,50,60,70,100]:

            #generate new workspace
            generateWorkspace.main(['-c',"cat==121 && l1pt>{0} && l2pt>{1} && bosonpt>{2}".format(l1pt,l2pt,bosonpt)])

            #run limits
            os.system('text2workspace.py shapes-parametric-140murad.datacard.dat')
            os.system('combine shapes-parametric-140murad.datacard.dat.root -M AsymptoticLimits -t -1')

            #read 50% quantile 
            fIn=ROOT.TFile.Open('higgsCombineTest.AsymptoticLimits.mH120.root')
            limit=fIn.Get('limit')
            limit.GetEntry(2)
            table.append((l1pt,l2pt,bosonpt,limit.limit))
            fIn.Close()


with open('optimresults.dat','w') as cache:
    for l in table:
        cache.write('%d & %d & %d & %3.2f\n'%l)
