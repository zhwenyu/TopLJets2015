#!/usr/bin/env python
"""
    launch a task separately per crossing angle
    add #raw muons to the ROC curves
"""

import numpy as np
import optparse
import os
import sys
import pickle
from matplotlib import pyplot as plt

def plotFeatureImportance(clf,features,plt):
    
    """ plots the feature importance """
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    y_pos=np.arange(len(features))
    ax.barh(y_pos,clf.feature_importances_,align='center',alpha=0.5)
    plt.yticks(y_pos,features)

    plt.xlabel('Feature importance')
    fig.text(0.12, 0.9, r'CMS', fontweight='bold')
    fig.text(0.18, 0.9, r'Preliminary')
    plt.show()

def showConfusionMatrix(clf,features,data,truth='class'):
    
    """displays the confusion matrix"""
    
    preds = clf.predict(data[features])
    pd.crosstab(data[truth], preds, rownames=['Observed'], colnames=['Predicted'])

def fitModels(data,features,opt):

    """ fits different models to the data """

    import pandas as pd
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
    from sklearn.model_selection import GridSearchCV
    import pickle

    best_models={'lin':{},'rfc':{},'gbc':{},'dnn':{}}
    
    xangle_list=data['s'][:,0]
    for xangle in [120.,130.,140.,150.]:

        filt=(xangle_list==xangle)
        X=data['X'][filt]
        y=data['y'][filt]
        print '[fitModels] @ crossing angle={0} murad has {1} events'.format(xangle,len(y))
        
        #prepare the train (and test) dataset    
        df=pd.DataFrame(X)
        df.columns=features
        df['class']=y.tolist()
        train=df.sample(frac=opt.trainFrac,random_state=200)
        test=df.drop(train.index)
        print train
        continue

        #optimize a random forest classifier
        rfc = RandomForestClassifier()
        rfc_params = {
            'bootstrap': [True],
            'n_estimators': [50,100], #,150,200],
            'max_depth': [2,5,], #10,15],
            'max_features':['sqrt'] #,'log2']
            }
        rfc_optim = GridSearchCV(estimator = rfc, param_grid = rfc_params, cv = 3, n_jobs = -1, verbose = 2)
        rfc_optim.fit(train[features], train['class'])
        print rfc_optim.best_estimator_.feature_importances_
        #bestFeatures=sorted(rfc_optim.best_estimator_.feature_importances_, key=lambda x: x[1])        
        best_models['rfc'][xangle]=(rfc_optim.best_estimator_,rfc_optim.best_params_,features)

        #optimize a gradient boost classifier
        gbc = GradientBoostingClassifier()
        gbc_params = {
            'loss' : ['deviance'],
            'learning_rate':[0.1], #0.01],
            'n_estimators':[50], #,100,200],
            'max_depth': [2], #,3,5],
            'max_features':['sqrt'], #,'log2']
            }
        gbc_optim = GridSearchCV(estimator = gbc, param_grid = gbc_params, cv = 3, n_jobs = -1, verbose = 2)
        gbc_optim.fit(train[branches], train['class'])
        best_models['gbc'][xangle]=(gbc_optim.best_estimator_,gbc_optim.best_params_)

    #save results  to pickle file
    out_url=os.path.join(opt.output,'pu_models.pck')
    with open(out_url,'w') as cache:    
        pickle.dump(best_models, cache, pickle.HIGHEST_PROTOCOL)
    print 'Best fit models have been stored in',out_url

def runTrainJob(url,features,spectators,categs,onlyThis,opt):


    """ 
    converts the ROOT trees to numpy arrays for the 
    - features to use in the training
    - the spectator variables 
    - the classes to predict
    a selection string is used to filter the original events in the trees
    """

    #build the chain
    from ROOT import TChain
    
    t=TChain('tree')
    for f in os.listdir(url):
        if onlyThis:
            if f!=onlyThis: continue
        else:
            if not 'DoubleMuon' in f and not 'SingleMuon' in f: 
                continue
        t.AddFile(os.path.join(url,f))
    print 'Data chain has {0} events'.format(t.GetEntries())

    #convert to numpy arrays
    from sklearn import preprocessing
    from root_numpy import tree2array

    data={}
    data['X']=tree2array(t, branches=features,   selection=opt.selection)
    data['s']=tree2array(t, branches=spectators, selection=opt.selection)
    data['y']=tree2array(t, branches=categs,     selection=opt.selection)

    #convert features and spectators to array of floats
    for tag,branches in [('X',features),('s',spectators)]:
        data[tag]=data[tag].astype([(b,'<f4') for b in branches])
        data[tag]=data[tag].view(np.float32).reshape(data[tag].shape + (-1,))

    #preprocess features
    data['X']=preprocessing.scale(data['X'])

    print 'Converted to numpy array'

    if opt.model is None:
        out_url=os.path.join(opt.output,'train_data.pck')
        with open(out_url,'w') as cache:    
            pickle.dump(data, cache, pickle.HIGHEST_PROTOCOL)
        print 'A dump of the data in numpy format is saved @',out_url
        fitModels(data,features,opt)
        


def runTrainJobPacked(args):
    
    """wrapper for parallel threads"""

    try:
        runTrainJob(*args)
    except Exception as e:
        print '>'*50
        print e
        print '>'*50
        return False

def main():


    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',
                      dest='input',   
                      default='/eos/cms/store/cmst3/user/psilva/ExclusiveAna/ab05162/Chunks/',
                      help='input directory with the files [default: %default]')
    parser.add_option('--trainFrac',
                      dest='trainFrac',
                      default=0.5,
                      type=float,
                      help='fraction to use for training [default: %default]')
    parser.add_option('-o', '--output',
                      dest='output',
                      default='test/analysis/pps/',
                      help='output directory [default: %default]')
    parser.add_option('-s', '--selection',
                      dest='selection',   
                      default='isZ && bosonpt<10 && (nRPtk==0 || nRPtk==2)',
                      help='selection [default: %default]')
    parser.add_option('-m', '--model',
                      dest='model',   
                      default=None,
                      help='pickle file with trained models [default: %default]')
    (opt, args) = parser.parse_args()

    #build the feature,spectator, category and selection
    features=['nvtx','rho']
    for v in ['PFMult','PFPz','PFHt','PFChMult','PFChPz','PFChHt']:
        for r in ['EB','EE','HE','HF']:
            if 'Ch' in v and 'H' in r : 
                continue
            for o in ['Sum','Diff']:
                features.append(v+o+r)
    spectators=['beamXangle','nrawmu-2']  #keep beamXangle in the first position...
    categs='nRPtk>0? 0. : 1.'

    #fit directly the data if pickle file is available
    if '.pck' in opt.input:
        with open(opt.input,'r') as cache:
            data=pickle.load(cache)
            print data['y']
        fitModels(data,features,opt)


    #build the tasks
    task_list=[]
    buildTrainData=True
    if opt.model is None:
        task_list.append( (opt.input, features, spectators, categs, None, opt) )
    else:
        task_list=[ (opt.input, features, spectators, categs, f, opt) for f in os.listdir(opt.input) ]

    #run it
    import multiprocessing as MP
    pool = MP.Pool(8)
    pool.map(runTrainJobPacked, task_list)

if __name__ == "__main__":
    sys.exit(main())


"""
# In[2]:


#import tree from file


print 'Preselection for training is',selection

inF = ROOT.TFile('/eos/user/p/psilva/Data13TeV_2017F_DoubleMuon.root','READ')
tree = inF.Get('tree')

#convert to all floats, apply standard scaling, and reshape as image
X=rnp.tree2array(tree,
                 branches=branches,
                 selection=selection)
X=X.astype([(b,'<f4') for b in branches])
X=X.view(np.float32).reshape(X.shape + (-1,))
X=preprocessing.scale(X)

#categories
y=rnp.tree2array(tree,
                 branches="nRPtk==0? 1. : 0.",
                 selection=selection)

print len(X),'events read'
print len(y[y>0]),'signal-like',len(y[y==0]),'background-like'

#all done
inF.Close()


# In[3]:


# In[36]:


#show correlation
corr_bkg=df.loc[df['class']==0].corr(method='pearson')
corr_sig=df.loc[df['class']==1].corr(method='pearson')
corr_ratio=corr_bkg/corr_sig

fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(corr_ratio,cmap=plt.cm.jet,vmin=-2,vmax=2)
ax.xaxis.tick_bottom()
clb=fig.colorbar(cax)
clb.set_label(r'$\rho$(B) / $\rho$(S)',rotation=270.,labelpad=15)
plt.xlabel('Variable index')
plt.ylabel('Variable index')
fig.text(0.23, 0.9, r'CMS', fontweight='bold')
fig.text(0.29, 0.9, r'Preliminary')
plt.show()


# In[6]:




# In[23]:


from keras.models import Sequential
from keras.layers.core import Dense, Dropout
from keras.optimizers import Adam

mlp = Sequential()
mlp.add(Dense(256, input_dim=len(branches), activation='relu'))
mlp.add(Dropout(0.1))
mlp.add(Dense(64, activation='relu'))
mlp.add(Dropout(0.1))
mlp.add(Dense(16, activation='relu'))
mlp.add(Dropout(0.1))
mlp.add(Dense(1, activation='sigmoid'))

mlp.compile(loss='binary_crossentropy',
              optimizer='rmsprop',
              metrics=['accuracy'])

mlp.fit(train[branches], train['class'],
          epochs=100,
          batch_size=128)


# In[7]:


    
    
with open('bdt_optim.pck','r') as cache:
    train_results=pickle.load(cache)

for key in train_results:
    print key,'summary'
    print train_results[key][1]
    plotFeatureImportance(clf=train_results[key][0],features=branches,plt=plt)
    showConfusionMatrix(clf=train_results[key][0],features=branches,data=test)


# In[30]:


#final summary
from sklearn.metrics import auc, roc_curve

fig=plt.figure(111)
plt.plot([0,1], [1,0], 'k--', color='orange')

for key in train_results:
    clf=train_results[key][0]
    y_prob=clf.predict_proba(test[branches])[:,0]
    fpr, tpr, thresholds = roc_curve(test['class'],y_prob, pos_label=0)
    aucVal = auc(fpr, tpr)
    plt.plot(1-fpr, tpr, label='{0} auc = {1:.3f}'.format(key,aucVal))


for v in ['nvtx','rho']:
    fpr,tpr,thresholds=roc_curve(test['class'],test[v],pos_label=0)
    aucVal=auc(fpr,tpr)
    plt.plot(1-fpr, tpr, label='{0} auc = {1:.3f}'.format(v,aucVal))
    

plt.xlabel('Signal efficiency')
plt.ylabel('Background rejection')
plt.legend(loc='lower left')
plt.grid()
fig.text(0.12, 0.9, r'CMS', fontweight='bold')
fig.text(0.18, 0.9, r'Preliminary')
plt.show()


# In[ ]:




"""
