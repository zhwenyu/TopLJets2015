#!/usr/bin/env python

import ROOT
import sys

fIn=ROOT.TFile.Open(sys.argv[1])
tables={'EM':[ (1,'visited'),
               (3,'dileptons'),
               (6,'>=2j'),
               (9,'>=2b') ],
        'EE':[ (1,'visited'),
               (3,'dileptons'),
               (4,'Z veto'),
               (5,'MET>40'),
               (6,'>=2j'),
               (9,'>=2b') ]}
tables['MM']=tables['EE']

for ch in tables:
    h=fIn.Get('cutflow_%s'%ch)
    print 'Event selection for',ch
    for bin,cut in tables[ch]:
        print '%15s'%cut,'%15d'%h.GetBinContent(bin)


fIn.Close()
