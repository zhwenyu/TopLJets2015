#! /usr/bin/env python
import os, sys

try:
    inputdir = sys.argv[1]
    outputdir = sys.argv[2]
    if not os.path.isdir(inputdir):
        print("Input directory not found:", inputdir)
        exit(-1)
except IndexError:
    print("Usage: mergeProduction.py inputdir outputdir")
    exit(-1)

for subdir in [x[0] for x in os.walk(inputdir)]:
    if (subdir == inputdir): continue

    outputsubdir = outputdir + '/' + subdir.rsplit('/', 1)[1]
    os.system('mkdir -p %s'%(outputsubdir))
    
    files = []
    for f in [x[2] for x in os.walk(subdir)][0]:
        files.append(subdir+'/'+f)
    chunks = [files[i:i + 10] for i in xrange(0, len(files), 10)]
    
    counter = 0
    for chunk in chunks:
        cmd = 'hadd %s/MergedMiniEvents_%i.root %s'%(outputsubdir, counter, ' '.join(chunk))
        print(cmd)
        os.system(cmd)
        counter += 1
