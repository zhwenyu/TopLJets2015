#!/usr/bin/env python

import pickle
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import pandas as pd
from keras.models import load_model
from matplotlib import pyplot as plt

def plotFeatureImportance(models,outName):
    
    """ plots the feature importance """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    features=[]
    for _,m_feats in models:
        for f,_ in m_feats:
            features.append(f)
    features=list(set(features))

    feature_importance={}
    for m,m_feats in models:
        feature_importance[m]=[0.]*len(features)
        for f,fval in m_feats:
            idx=features.index(f)
            feature_importance[m][idx]=fval

    y_pos=np.arange(len(features))
    for m in feature_importance:
        ax.barh(y_pos,feature_importance[m],align='center',alpha=0.5,label=m)
    plt.yticks(y_pos,features)
    plt.tick_params(axis='both', which='major', labelsize=8)
    plt.tick_params(axis='both', which='minor', labelsize=8)

    plt.xlabel('Feature importance')
    fig.text(0.20, 0.96, r'CMS', fontweight='bold')
    fig.text(0.26, 0.96, r'Preliminary')
    plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.1)
    plt.legend(loc='upper right')
    plt.grid()
    plt.savefig(outName)    
    plt.clf()

def showConfusionMatrix(clfList,data,truth,outfile):
    
    """displays the confusion matrix"""

    from sklearn.metrics import confusion_matrix

    for m,clf,features in clfList:
        if isinstance(clf,basestring) : clf=load_model(clf)     
        y_pred = clf.predict(data[features])
        conf = confusion_matrix(data[truth], y_pred)


def showROCs(clfList,data,truth,outfile):

    """ show final ROC curves"""

    from sklearn.metrics import auc, roc_curve

    fig=plt.figure(111)
    plt.plot([0,1], [0,1], ':', color='gray')
    plt.plot([0.25,0.25], [0,0.4], linewidth=2, color='gray')

    for m,clf,features in clfList:
        try:
            y_prob = clf.predict_proba(data[features])[:,0]
        except:
            continue
            pass
        fpr, tpr, thresholds = roc_curve(data[truth],y_prob,pos_label=0)
        aucVal = 1-auc(tpr, fpr)
        plt.plot(tpr, fpr, label='{0} auc = {1:.3f}'.format(m,aucVal))


    for v in ['nvtx','rho','nrawmu']:
        fpr,tpr,thresholds=roc_curve(data['class'],data[v],pos_label=0)
        fpr=1-fpr # these variables are cut by upper value v<max_val => signal-like
        tpr=1-tpr
        aucVal=1-auc(tpr,fpr)
        plt.plot(tpr, fpr, linestyle='--', label='{0} auc = {1:.3f}'.format(v,aucVal))
    
    plt.xlabel('Signal efficiency')
    plt.ylabel('Background efficiency')
    plt.legend(loc='upper left')
    plt.grid()
    fig.text(0.12, 0.9, r'CMS', fontweight='bold')
    fig.text(0.18, 0.9, r'Preliminary')
    plt.savefig(outfile)    
    plt.clf()

with open('test/analysis/pps/pu_models.pck','r') as cache:
    best_models=pickle.load(cache)

with open('test/analysis/pps/train_data.pck','r') as cache:
    data=pickle.load(cache)
    features=pickle.load(cache)
    spectators=pickle.load(cache)


xangle_list=data['s'][:,0]
for xangle in [120,130,140,150]:

    clfList=[]
    rfcList=[]
    for key in best_models:
        clfList.append( (key,best_models[key][xangle][0],best_models[key][xangle][-1]) )
        if not 'rfc' in key: continue
        rfcList.append( (key,best_models[key][xangle][2]) )

    #feature importance (for random forest classifier only)
    plotFeatureImportance(rfcList,'featimportance_%d.png'%xangle)
    
    #filter the data for this crossing angle
    filt=(xangle_list==xangle)
    X=data['X'][filt]
    s=data['s'][filt][:,-1]
    y=data['y'][filt]
    df=pd.DataFrame(X)
    df.columns=features
    df['nrawmu']=s.tolist()
    df['class']=y.tolist()

    #showConfusionMatrix(clfList=clfList,data=df,truth='class',outfile='confmatrix_%d.png'%xangle)
    showROCs(clfList=clfList,data=df,truth='class',outfile='rocs_%d.png'%xangle)

