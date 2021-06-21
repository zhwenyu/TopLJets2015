import pickle
import ROOT

with open('TT2l_npsysts.pck','rb') as fin:
    np_weights=pickle.load(fin)

cat=('2b','hpt')
var_list=['CRerd','CRgmove','CRqcd','UEup','UEdown','hdampup','hdampdown']

print('Weights for cat=',cat)
print('m(l,b)\t {}'.format('\t\t'.join(var_list)))
print('='*110)
for mlb in [20,60,100,150,200]:
    line='{:d}\t'.format(mlb)
    for var in var_list:
        line += '{:3.2f}/{:3.2f}\t'.format(np_weights[var][cat]['Up'].Eval(mlb),np_weights[var][cat]['Down'].Eval(mlb))
    print line
