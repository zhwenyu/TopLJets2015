import ROOT
import os
import generateBinnedWorkspace
import sys
import pickle

m=int(sys.argv[1])
xangle=int(sys.argv[2])
baseDir=sys.argv[3]

#check which points are available
with open('%s/scanpoints.pck'%baseDir,'r') as cache:
    task_list=pickle.load(cache)

#scan results
hr95=ROOT.TH1F('optim_r95',';Scan point;95% CL limit on #sigma/#sigma_{th}',len(task_list),0,len(task_list))
hsig=hr95.Clone('optim_sig')
hsig.GetYaxis().SetTitle('Expected significance (asymptotic)')
scanResults=[]

#run combine at each point
nTasks=len(task_list)
for task in task_list:

    itask=task[0]
    print 'Starting',itask,'/',nTasks
    taskDir=baseDir+'/optim_task%d'%itask
    iwpResults=[]

    try:

        #run combine
        dc=[os.path.join(taskDir,f) for f in os.listdir(taskDir) if 'dat' in f and 'a%d'%xangle in f]
        combStr=' '.join( ['cat%d=%s'%(i,dc[i]) for i in range(len(dc))] )    
        os.system('combineCards.py {0} > combined_card.dat\n'.format(combStr))
        os.system('text2workspace.py combined_card.dat -o workspace.root\n')    
    
        #run limits and read 50% quantile
        os.system('combine workspace.root -M AsymptoticLimits -t -1\n')
        fIn=ROOT.TFile.Open('higgsCombineTest.AsymptoticLimits.mH120.root')
        limit=fIn.Get('limit')
        for ientry in range(5):
            limit.GetEntry(ientry)
            iwpResults.append(limit.limit)
        r95=iwpResults[2]
        hr95.SetBinContent(itask,r95)
        fIn.Close()

        #run significance
        os.system('combine workspace.root -M Significance -t -1 --expectSignal=1\n')
        fIn=ROOT.TFile.Open('higgsCombineTest.Significance.mH120.root')
        limit=fIn.Get('limit')
        limit.GetEntry(0)
        exp_sig=limit.limit
        iwpResults.append(exp_sig)
        hsig.SetBinContent(itask,exp_sig)
        fIn.Close()

        scanResults.append( [x for x in task]+[r95,exp_sig] )
    except:
        print 'Failed to analyse results for',task    

with open('{0}/optim_results_m{1}_a{2}.pck'.format(baseDir,m,xangle),'w') as cache:
    pickle.dump(scanResults,cache,pickle.HIGHEST_PROTOCOL) 
    pickle.dump(hr95,cache,pickle.HIGHEST_PROTOCOL) 
    pickle.dump(hsig,cache,pickle.HIGHEST_PROTOCOL) 
