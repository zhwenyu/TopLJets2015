import os

eras=[('B',297020,299329),
      ('C1',299337,300785),
      ('C2',300806,302029),
      ('D',302030,303434),
      ('E',303435,304826),
      ('F1',304911,305114),
      ('F2',305178,305902),
      ('F3',305965,306462),]

for era,rmin,rmax in eras:
    print(era,rmin,rmax)
    os.system('filterJSON.py  Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt --min {} --max {} --output Cert_{}.json'.format(rmin,rmax,era))
    os.system('brilcalc lumi -b "STABLE BEAMS" -u /pb -i Cert_{0}.json --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json  > lumi_{0}.dat'.format(era))
