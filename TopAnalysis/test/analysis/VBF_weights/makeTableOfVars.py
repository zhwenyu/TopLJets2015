import os,sys
from subprocess import Popen, PIPE
import re

VARDICT={
"centralEta"  : ('$\\eta_{cen}$','most central $\\eta$'),
"forwardeta"  : ('$\\eta_{fwd}$','most forward $\\eta$'),
"leadj_pt"    : ('$p_{T,1}$','leading jet $p_T$'),
"subleadj_pt" : ('$p_{T,2}$','sub-leading jet $p_T$'),
"mjj"        : ('$m_{jj}$','dijet invariant mass'),
"detajj"     : ('$\\delta\\eta_{jj}$','dijet rapidity opening'),
"jjpt"       : ('$p_{T,jj}$','dijet $p_T$'),
"dphijj"     : ('$\\delta\\phi_{jj}$','dijet azimuthal angle opening'),
"ystar"      : ('$y^*$','boson rapidity in dijet rest frame'),
"relbpt"     : ('$p_{T,rel}$','relative boson $p_T$ with respect to vis. system'),
"dphibjj"      : ('$\\Delta\phi(b,jj)$','azimuthal angle between boson and dijet'),
"balance"      : ('$p_{T,bal}$','balance of the vis. system'),
"leadj_gawidth"   : ('$w_1$','lead jet width'),
"subleadj_gawidth": ('$w_2$','sub-lead jet width'),
"subleadj_c2_02"  : ('$C_2(2,2)$','energy correlation function for sub-lead jet'),
"jjetas"          : ('$\\eta_{jj}$','rapidity of the dijet system'),
"centjy"          : ('$y_3$','rapidity of the third jet'),
"ncentjj"         : ('$N_{j,cen}$','multiplicity of jets in rapidity system'),
"dphivj0"         : ('$\\Delta\\phi(b,j_1)$', 'azimuthal angle between boson and jet 1'),
"dphivj1"         : ('$\\Delta\\phi(b,j_2)$', 'azimuthal angle between boson and jet 2'),
"dphivj2"         : ('$\\Delta\\phi(b,j_3)$', 'azimuthal angle between boson and jet 3'),
"dphivj3"         : ('$\\Delta\\phi(b,j_4)$', 'azimuthal angle between boson and jet 4'),
"mht"             : ('MHT','missing $H_T$'),
"ht"              : ('$H_T$','scalar $H_T$'),
"isotropy"        : ('isotropy','event isotropy'),
"circularity"     : ('circularity','event circularity'),
"sphericity"      : ('sphericity','event sphericity, $1.5*(q_1+q_2)$'),
"aplanarity"      : ('aplanarity','event aplanarity, $1.5*q_1$'),
"C"               : ('C','event C, $3.*(q_1*q_2+q_1*q_3+q_2*q_3)$'),
"D"               : ('D','event D, $27.(q_1.q_2.q_3)$'),
"jet_c2_001"	  : ('$C_1(2,2)$','energy correlation function for leading jet'),
"jet_c2_002"	  : ('$C_2(2,2)$','energy correlation function for sub-lead jet'),
"subleadjet_qg"	  : ('$d_{qg}^2$','quark-gluon discriminator of sub-lead jet'),
"leadjet_qg"	  : ('$d_{qg}^1$','quark-gluon discriminator of leading jet'),
}


data=Popen(["grep -ir Variable %s | awk '{print $5}'"%sys.argv[1]], stdout=PIPE, shell=True)
out,err=data.communicate()

print '\\begin{table}[!ht]'
print '\\topcaption{Variables as input to the MVA discriminator. In event shape variables, $q_i$s are eigenvalues of the momentum tensor $\\frac{\sum{p_j[a]*p_j[b]}}{\sum{p_j^2}}$ normalized to one.} '
print '\centering'
print '\\begin{tabular}{|c|c|} \hline'
pattern = r'"([A-Za-z0-9_\./\\-]*)"'
for line in out.split():
    if not 'Title' in line : continue
    m = re.search(pattern, line).group().replace("\"","")
    if m in VARDICT:
        print '%25s & %50s \\\\'%VARDICT[m]
    else:
        print '%25s & %50s \\\\'%(m,'')
print '\hline'
print '\end{tabular} '
print '\end{table}'
