import os
import ROOT
import sys
import optparse
from TopLJets2015.TopAnalysis.storeTools import getEOSlslist

"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',       dest='inDir',       help='input directory with files',               default=None,   type='string')
    parser.add_option('-o', '--outDir',      dest='outDir',      help='output directory with files',              default=None,   type='string')
    parser.add_option('-c', '--cleanup',     dest='cleanup',     help='removes original crab directory',          default =False, action='store_true')
    parser.add_option(      '--nocheck',     dest='nocheck',     help='do not prompt user',                       default=False,  action='store_true')
    parser.add_option(      '--only',        dest='only',        help='only this tag',                            default=None,   type='string')
    (opt, args) = parser.parse_args()

    baseEOS='/eos/cms/'

    #prepare output directory
    if opt.outDir is None: opt.outDir=opt.inDir
    os.system('mkdir -p %s/%s'%(baseEOS,opt.outDir))

    dset_list=getEOSlslist(directory=opt.inDir,prepend='')
    for dset in dset_list:
        dsetname=dset.split('/')[-1]

        pub_list=getEOSlslist(directory=dset,prepend='')
        for pubDir in pub_list:

            if not 'crab' in pubDir:
                print 'Ambiguity found @ <publication-name> for <primary-dataset>=%s , bailing out'%dsetname
                continue

            #select if it doesn't match the required selection
            pub=pubDir.split('/crab_')[-1]
            if opt.only:
                if pub!=opt.only: 
                    continue

            #check if it's an extension
            pubExt=None
            try:
                extSplit=pub.split('_ext')
                pubExt='ext%d'%(len(extSplit)-1)
                pub=extSplit[0]
                print 'Extension will be postfixed with ',pubExt
            except:
                print 'Core sample (no extension)'
                
            time_list=getEOSlslist(directory=pubDir,prepend='')
            if len(time_list)!=1:
                print 'Ambiguity found @ <time-stamp> for <primary-dataset>=%s , bailing out'%dsetname
                continue
            time_stamp=time_list[0].split('/')[-1]

            out_list=[]
            count_list=getEOSlslist(directory=time_list[0],prepend='')
            for count in count_list: 
                if '/merge' in count and  'group/hintt' in count : continue
                out_list += getEOSlslist(directory=count,prepend='')
            file_list=[x for x in out_list if '.root' in x]

            newDir='%s/%s' % (opt.outDir,pub)        
            print '<primary-dataset>=%s <publication-name>=crab_%s <time-stamp>=%s has %d files' % (dsetname,pub,time_stamp,len(file_list) )
            if not opt.nocheck:
                choice = raw_input('Will move to %s current output directory. [y/n] ?' % newDir ).lower()
                if not 'y' in choice : continue
            
            os.system('mkdir -p %s/%s'%(baseEOS,newDir))
    
            moveIndividualFiles=True
            if len(file_list)>0:
                subgroupMerge=10 
                if 'Data' in pub: subgroupMerge=50
                moveIndividualFiles=False

                splitFunc = lambda A, n=subgroupMerge: [A[i:i+n] for i in range(0, len(A), n)]
                split_file_lists = splitFunc( file_list )
                
                for ilist in xrange(0,len(split_file_lists)):
                    if pubExt:
                        mergedFileName='MergedMiniEvents_%d_%s.root '%(ilist,pubExt)
                    else:
                        mergedFileName='MergedMiniEvents_%d.root '%ilist
                    toAdd='%s ' % mergedFileName                    
                    for f in split_file_lists[ilist]:                            
                        toAdd += '%s/%s '%(baseEOS,f) 
                    finalOutput='%s/%s/%s'%(baseEOS,newDir,mergedFileName)
                    fIn=ROOT.TFile.Open(finalOutput)
                    try:
                        if fIn or not fIn.IsZombie():
                            print '%s already in EOS, skipping'%finalOutput
                            fIn.Close()
                            continue
                    except:
                        pass
                        
                    os.system('hadd -f -k %s'%toAdd)
                    os.system('mv -v %s %s/%s/%s'%(mergedFileName,baseEOS,newDir,mergedFileName))
                        
            #if still needed copy individual files
            if moveIndividualFiles:
                for f in file_list : 
                    newF=f
                    if pubExt:
                        newF=f.replace('.root','_%s.root'%pubExt)
                    os.system('mv -v %s %s/%s/%s'%(f,baseEOS,newDir,newF))

            if not opt.nocheck and opt.cleanup : 
                choice = raw_input('Will remove output directory. [y/n] ?').lower()
                if 'y' in choice: 
                    os.system('rm -r %s/%s'%(baseEOS,dset))

            print 'Crab outputs may now be found in %s' % newDir

    print '-'*50
    print 'All done. In case errors were found check that the crab output structure is '
    print '<outLFNDirBase>/<primary-dataset>/<publication-name>/<time-stamp>/<counter>/<file-name>'
    print '-'*50
        


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

