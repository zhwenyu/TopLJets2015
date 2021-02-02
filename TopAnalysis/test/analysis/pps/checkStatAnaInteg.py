import ROOT
import sys
import os

from prepareOptimScanCards import MASSPOINTS,OPTIMLIST

def checkTemplates(fIn,b,massList):

    missMasses=[]
    for m in massList:
        try:
            h=fIn.Get(('fidsig_{0}_m{1}'.format(b,m)))
            h.Integral()
        except Exception as e:
            missMasses.append(m)

    return missMasses



_forceRecreateDataCards=True

baseDir=sys.argv[1]
optimList=[os.path.join(baseDir,d) for d in os.listdir(baseDir) if 'optim_' in d] 
toCheck=[]
for d in optimList:

    optimIdx=int(d.split('_')[-1])

    for ch in [121,169,22]:

        allGoodToCombine=True
        for tag in ['bkg']+['_'.join( [str(y) for y in x] ) for x in MASSPOINTS]:

            try:
                fIn=ROOT.TFile.Open(os.path.join(d,'shapes_%d_%s.root'%(ch,tag)))
                if fIn.IsZombie():
                    raise Exception('corrupted')
                elif fIn.TestBit(ROOT.TFile.kRecovered) : 
                    raise Exception('recovered')
                else:
                    size=fIn.GetListOfKeys().GetSize()
                    if size==0: 
                        raise Exception('no keys')

                #check histograms are there
                b='zee'
                if ch==169:  b='zmm'
                if ch==22:   b='g'
                if tag!='bkg':
                    
                    missMasses=checkTemplates(fIn,b,tag.split('_'))
                    if len(missMasses)>0:
                        
                        newTag='_'.join(missMasses)
                        newUrl=os.path.join(d,'shapes_%d_%s.root'%(ch,newTag))
                        if os.path.isfile(newUrl):
                            
                            print '*'*50
                            print 'Trying to find',missMasses,'in recovery job',newUrl
                            print '*'*50

                            newFin=ROOT.TFile.Open(newUrl)
                            missMasses=checkTemplates(newFin,b,missMasses)
                            newFin.Close()
                        
                        if len(missMasses)>0:
                            allGoodToCombine=False
                            print 'Will submit',d,'for',ch,'tag=',tag,'miss massess:',missMasses
                            toCheck.append( '%d,%d,0,\\"%s\\"'%(optimIdx,ch,','.join(missMasses)) )
                else:
                    try:
                        h=fIn.Get("bkg_"+b)
                        h.Integral()
                    except Exception as e:
                        print '*'*50
                        print tag,e
                        print '*'*50
                        toCheck.append( '%d,%d,1,\\"-1\\"'%(optimIdx,ch) )

                fIn.Close()

            except Exception as e:
                
                #corrupted file? submit all and remove file

                print('Removing corrupted file')
                os.system('rm -v {}'.format(os.path.join(d,'shapes_%d_%s.root'%(ch,tag))))

                allGoodToCombine=False
                print 'Will submit',d,'for',ch,'tag=',tag,'error:',e
                doBackground=1 
                massList=-1
                if tag!='bkg': 
                    doBackground = 0
                    massList     = tag.replace('_',',')
                toCheck.append( '%d,%d,%d,\\"%s\\"'%(optimIdx,ch,doBackground,massList) )


        if not allGoodToCombine:
            print 'Some masses missing ({}) but will still combine'.format(missMasses)
            #continue

        if os.path.isfile(os.path.join(d,'shapes_{}.root'.format(ch))) and not _forceRecreateDataCards:
            print ch,'already combined for',d
            continue

        print 'Will combine',d,'for',ch
        cuts = OPTIMLIST[ optimIdx ]
        os.system('python test/analysis/pps/generateBinnedWorkspace.py --doDataCards -o {} --finalStates {} --cuts "{}"'.format(d,ch,cuts))

if len(toCheck)>0:

    resub=os.path.join(baseDir,'zxstatana_scan_resub.sub')
    with open(resub,'w') as f:
        f.write('executable   = %s/optim_$(optimId)/optimJob_$(finalState).sh\n'%os.path.abspath(baseDir))
        f.write('arguments    = $(doBackground) $(massList)\n')
        f.write('output       = zxstatana_scan_resub.out\n')
        f.write('error        = zxstatana_scan_resub.err\n')
        f.write('log          = zxstatana_scan_resub.log\n')
        f.write('+JobFlavour  = "tomorrow"\n')
        f.write('request_cpus = 4\n')
        f.write('queue optimId,finalState,doBackground,massList from (\n')
        for j in toCheck: f.write('\t%s\n'%j)
        f.write(')\n')

    print '%d resubmission jobs can be found in %s'%(len(toCheck),resub)
    #os.system('condor_submit %s'%resub)
