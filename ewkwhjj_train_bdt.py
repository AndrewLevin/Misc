#!/usr/bin/env python
# coding: utf-8

# In[1]:

import xgboost as xgb
from coffea.util import load
import coffea.hist
import pandas
import matplotlib.pyplot as plt
import mplhep as hep
from sklearn.model_selection import train_test_split 
import numpy as np
from coffea import hist
import time
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--gpu',dest='gpu',type=int,default=0)

args = parser.parse_args()

# In[2]:


#result_2016_pre = load("training_data/training_data_2016pre")
#result_2016_post = load("training_data/training_data_2016post")
#result_2017 = load("training_data/training_data_2017")
#result_2018 = load("training_data/training_data_2018")
result_2016pre = load("training_data_2016pre")
result_2016post = load("training_data_2016post")
result_2017 = load("training_data_2017")
result_2018 = load("training_data_2018")
#result_2018 = load("/afs/cern.ch/project/afs/var/ABS/recover/R.1935065265.07291343/tmp/training_data_2018")

resolved_variables = ['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','met','metphi','higgsjet1pt','higgsjet2pt','vbfjet1pt','vbfjet2pt','higgsjet1eta','higgsjet2eta','vbfjet1eta','vbfjet2eta','higgsjet1phi','higgsjet2phi','vbfjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbfjet1btag','vbjfet2btag','higgsdijetmass','vbfdijetmass','vbfdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt','higgdijetdeltar','vbfdijetdeltar','higgsdijetabsdeta','dr1','dr2','dr3','dr4']

xs = {
    'qcdwphjj' : 0.1154*0.5840*0.33, 
    'qcdwmhjj' : 0.06917*0.5840*0.33,
    'qcdwph' : 0.2832*0.5840, 
    'qcdwmh' : 0.1770*0.5840,
    'ttsemi' : 365.4,
    'tthad' : 377.96,
    'ewkwhjj_reweighted' : 0.02656*0.5840,
    'ewkwhjj' : 0.02656*0.5840,
    'wlep' : 61526.7,
    'wlepsherpa' : 61526.7,
    'wlep2j' : 3338.0,
    'wlepfxfx0j' : 53330.0,
    'wlepfxfx1j' : 8875.0,
    'wlepfxfx2j' : 3338.0,
    'wlepmlm1j' : 8873.0,
    'wlepmlm2j' : 2793.0,
    'wlepmlm3j' : 992.5,
    'wlepmlm4j' : 380,
    'wlepht70100' : 1264.0,
    'wlepht100200' : 1256.0,
    'wlepht200400' : 335.5,
    'wlepht400600' : 45.25,
    'wlepht600800' : 10.97,
    'wlepht8001200' : 4.933,
    'wlepht12002500' : 1.16,
    'wlepht2500' : 0.008001,
    'wbb' : 215,
    'stoptchan' : 136.02,
    'santitoptchan' : 80.95
}

lumi = {'2016pre' : 19.5, '2016post' : 16.8, '2017' : 41.5, '2018' : 59.6}


# In[3]:


training_data_resolved_2016pre = pandas.concat(
    [                                                                                                                                                                                      
        pandas.concat(                                                                                                                                                                     
            [                                                                                                                                                                              
                pandas.DataFrame(result_2016pre["variables"][key].value,columns=resolved_variables),                                                               
                pandas.DataFrame(result_2016pre["weights"][key].value,columns=["weight"]),                                                                                                         
                pandas.DataFrame(len(result_2016pre["variables"][key].value)*[key],columns=["label"])                                                                                              
         ],axis=1)                                                                                                                                                                         
        for key in result_2016pre['variables'].keys()                                                                                                                                          
    ],                                                                                                                                                                                     
ignore_index=True)

training_data_resolved_2016post = pandas.concat(                                                                                                                                                                        
    [                                                                                                                                                                                      
        pandas.concat(                                                                                                                                                                     
            [                                                                                                                                                                              
                pandas.DataFrame(result_2016post["variables"][key].value,columns=resolved_variables),                                                               
                pandas.DataFrame(result_2016post["weights"][key].value,columns=["weight"]),                                                                                                         
                pandas.DataFrame(len(result_2016post["variables"][key].value)*[key],columns=["label"])                                                                                              
         ],axis=1)                                                                                                                                                                         
        for key in result_2016post['variables'].keys()                                                                                                                                          
    ],                                                                                                                                                                                     
ignore_index=True)

training_data_resolved_2017 = pandas.concat(                                                                                                                                                                        
    [                                                                                                                                                                                      
        pandas.concat(                                                                                                                                                                     
            [                                                                                                                                                                              
                pandas.DataFrame(result_2017["variables"][key].value,columns=resolved_variables),                                                               
                pandas.DataFrame(result_2017["weights"][key].value,columns=["weight"]),                                                                                                         
                pandas.DataFrame(len(result_2017["variables"][key].value)*[key],columns=["label"])                                                                                              
         ],axis=1)                                                                                                                                                                         
        for key in result_2017['variables'].keys()                                                                                                                                          
    ],                                                                                                                                                                                     
ignore_index=True)

training_data_resolved_2018 = pandas.concat(                                                                                                                                                                        
    [                                                                                                                                                                                      
        pandas.concat(                                                                                                                                                                     
            [                                                                                                                                                                              
                pandas.DataFrame(result_2018["variables"][key].value,columns=resolved_variables),                                                               
                pandas.DataFrame(result_2018["weights"][key].value,columns=["weight"]),                                                                                                         
                pandas.DataFrame(len(result_2018["variables"][key].value)*[key],columns=["label"])                                                                                              
         ],axis=1)                                                                                                                                                                         
        for key in result_2018['variables'].keys()                                                                                                                                          
    ],                                                                                                                                                                                     
ignore_index=True)


# In[4]:


def apply_xs(training_data,result,lumi):
        
    weights = []
    
    for i in range(len(training_data)):
        assert(training_data["label"][i] in xs.keys())    

        weights.append(xs[training_data["label"][i]]*1000*lumi/result["nevents"][training_data["label"][i]])
        
    weights=pandas.DataFrame(weights,columns=["weight"])
    
    training_data = training_data[training_data.columns.drop(["weight"])]

    training_data = pandas.concat([training_data,weights],axis=1)
    
    return training_data

training_data_resolved_2016pre=apply_xs(training_data_resolved_2016pre,result_2016pre,lumi['2016pre'])
training_data_resolved_2016post=apply_xs(training_data_resolved_2016post,result_2016post,lumi['2016post'])
training_data_resolved_2017=apply_xs(training_data_resolved_2017,result_2017,lumi['2017'])
training_data_resolved_2018=apply_xs(training_data_resolved_2018,result_2018,lumi['2018'])


# In[5]:



training_data = pandas.concat([training_data_resolved_2016pre,
               training_data_resolved_2016post,
               training_data_resolved_2017,
               training_data_resolved_2018])

#training_data = training_data_resolved_2018

training_data= training_data.reset_index(drop=True)

print('training_data',len(training_data))

print('ewkwhjj',len(training_data[training_data["label"] == "ewkwhjj"]))

print('ttsemi',len(training_data[training_data["label"] == "ttsemi"]))

print('wlep',len(training_data[training_data["label"] == "wlep"]))

print('singletop',len(training_data[training_data["label"].str.contains("chan")]))

print('wlepht',len(training_data[training_data["label"].str.contains("wlepht")]))

print('wlepmlm',len(training_data[training_data["label"].str.contains("wlepmlm")]))

print('wlepfxfx',len(training_data[training_data["label"].str.contains("wlepfxfx")]))

#training_data = training_data[(training_data["label"] == "ewkwhjj") | (training_data["label"] == "ttsemi")]
training_data = training_data[(training_data["label"] == "ewkwhjj")
                              | (training_data["label"] == "ttsemi")
                              | (training_data["label"] == "stoptchan")
                              | (training_data["label"] == "santitoptchan")
#                              | training_data["label"].str.contains("wlepht")
                              | (training_data["label"] == "wlepht200400")
                                 ]

training_data= training_data.reset_index(drop=True)

weights = []
labels = []

for i in range(len(training_data)):
    assert(training_data["label"][i] in xs.keys())    

    #if  "wlepht" in training_data["label"][i]:
    #    xs_weights.append(xs[training_data["label"][i]]*1000*lumi/(result_2016_post["nevents"][training_data["label"][i]]+result_2016_post["nevents"][training_data["label"][i]]+result_2017["nevents"][training_data["label"][i]]+result_2018["nevents"][training_data["label"][i]]))
    #else:
    weights.append(training_data['weight'][i])
        
    if training_data["label"][i] == "ewkwhjj":
        labels.append(1)
    elif training_data["label"][i] == "ttsemi":
        labels.append(0)
    else:    
        labels.append(0)
        
assert(len(weights) == len(labels))
assert(len(weights) == len(training_data))
    
weight_sums = {}   
    
for i in range(len(weights)):
    
    if labels[i] in weight_sums.keys():
        weight_sums[labels[i]] += weights[i]
    else:
        weight_sums[labels[i]] = weights[i]       
        
#for i in range(len(weights)):
#    weights[i] = 1000000000*weights[i]/weight_sums[labels[i]]   
    
##weight_sums[1] = weight_sums[1]*20
#weight_sums[1] = weight_sums[1]*10.
weight_sums[1] = weight_sums[1]
       
for i in range(len(weights)):
    weights[i] = 10000*weights[i]/weight_sums[labels[i]]     
    
#unique_weights = []        

#for i in range(len(weights)):
#    if weights[i] not in unique_weights:
#        unique_weights.append(weights[i])

#print('andrew debug 1')        
#print(unique_weights)
#print('andrew debug 2')    
    
    
labels=pandas.DataFrame(labels,columns=["label"])
weights=pandas.DataFrame(weights,columns=["weight"])
#weights=pandas.DataFrame(labels,columns=["weight"])
#weights=pandas.DataFrame(np.ones(len(weights)),columns=["weight"])

training_data = training_data[training_data.columns.drop(["weight","label"])]

training_data = pandas.concat([training_data,labels,weights],axis=1)


# In[6]:


X = training_data[training_data.columns.drop(["weight","label"])]

w = training_data["weight"]

#print(float(len(w[w == -1]))/len(w))

y = training_data["label"]

X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(X, y, w, test_size=0.1, random_state=1)
    
dtrain = xgb.DMatrix(X_train, label=y_train)

dtest = xgb.DMatrix(X_test, label=y_test)

from xgboost import XGBClassifier,XGBRegressor


#xbgr = XGBRegressor(n_estimators=400,learning_rate=0.1,verbosity=1, max_depth=6,objective='reg:squaredlogerror')
#xbgr = XGBRegressor(n_estimators=400,learning_rate=0.1,verbosity=1, max_depth=6)
#xbgr = XGBClassifier(n_estimators=400,learning_rate=0.1,verbosity=1, max_depth=6)

xbgr = XGBClassifier(n_estimators=100000,learning_rate=0.1,verbosity=1, max_depth=15, tree_method='gpu_hist', gpu_id=args.gpu)
#xbgr = XGBClassifier(n_estimators=100000,learning_rate=0.1,verbosity=1, max_depth=15)

##xbgr = XGBClassifier(n_estimators=100000,learning_rate=0.01,verbosity=1, max_depth=15, tree_method='gpu_hist')
##xbgr = XGBClassifier(n_estimators=400,learning_rate=0.1,verbosity=1, max_depth=8)

import os
os.system('date')
xbgr.fit(X_train,y_train,sample_weight=w_train,sample_weight_eval_set=[w_test],early_stopping_rounds=10,eval_set=[(X_test, y_test)])
os.system('date')
#%time xbgr.fit(X_train,y_train,early_stopping_rounds=10,eval_set=[(X_test, y_test)])
#xbgr.fit(X_train,y_train,sample_weight=w_train,sample_weight_eval_set=[w_test],early_stopping_rounds=10,eval_set=[(X_test, y_test)],eval_metric='logloss')
#xbgr.fit(X_train,y_train,sample_weight=w_train,sample_weight_eval_set=[w_test],early_stopping_rounds=10,eval_set=[(X_test, y_test)],eval_metric='logloss')
#xbgr.fit(X_train,y_train,early_stopping_rounds=10,eval_set=[(X_test, y_test)],eval_metric='logloss')
#xbgr.fit(X_train,y_train,early_stopping_rounds=10,eval_set=[(X_test, y_test)],eval_metric='auc')
#xbgr.fit(X_train,y_train,sample_weight=w_train,sample_weight_eval_set=[w_test],early_stopping_rounds=10,eval_set=[(X_test, y_test)])
#xbgr.fit(X_train,y_train,sample_weight=w_train,sample_weight_eval_set=[w_test],eval_set=[(X_test, y_test)])

xbgr.save_model("model.bin")



plt.style.use(hep.style.ROOT) 

importancefigax=xgb.plot_importance(xbgr)

importancefig=importancefigax.figure

importancefig.savefig('importance.png', bbox_inches='tight')

from sklearn.metrics import roc_curve, auc

print(xbgr.evals_result())

xbgr.predict(X_test,iteration_range=(0, 1))
xbgr.predict(X_test,iteration_range=(0, 2))
xbgr.predict(X_test,iteration_range=(0, 3))
xbgr.predict(X_test,iteration_range=(0, 4))
xbgr.predict(X_test,iteration_range=(0, 5))

#y_test_pred = xbgr.predict(X_test) 
y_test_pred = xbgr.predict_proba(X_test)[:,1]

test_false_positive_rates, test_true_positive_rates, _ = roc_curve(np.array(y_test[:]==1,dtype=int), y_test_pred[:])

print(auc(test_false_positive_rates, test_true_positive_rates))

#y_train_pred = xbgr.predict(X_train) 
y_train_pred = xbgr.predict_proba(X_train) [:,1]

train_false_positive_rates, train_true_positive_rates, _ = roc_curve(np.array(y_train[:]==1,dtype=int), y_train_pred[:])

print(auc(train_false_positive_rates, train_true_positive_rates))

plt.style.use(hep.style.ROOT) 

rocfig = plt.figure()

plt.plot(1-test_false_positive_rates, test_true_positive_rates)
plt.plot(1-train_false_positive_rates, train_true_positive_rates)

plt.xlabel('Background Rejection Rate')
plt.ylabel('Signal Efficiency')

rocfig.savefig('roc.png')

#for i in y_test.index:
#    w_test[i] = weight_sums[y_test[i]]*w_test[i]/1000000000.
    
for i in y_test.index:
    w_test[i] = weight_sums[y_test[i]]*w_test[i]/10000
    
for i in y_train.index:
    w_train[i] = weight_sums[y_train[i]]*w_train[i]/10000
        
#h=coffea.hist.Hist('Events',coffea.hist.Cat('dataset', 'Dataset'),coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1))
h=coffea.hist.Hist('Events',coffea.hist.Cat('dataset', 'Dataset'),coffea.hist.Bin('bdtscore', 'BDT score', 80, 0.9, 1))
h.fill(dataset='Signal test',bdtscore=y_test_pred[y_test == 1],weight=w_test[y_test==1])
h.fill(dataset='Background test',bdtscore=y_test_pred[y_test == 0],weight=w_test[y_test==0])
#h.fill(dataset='Signal train',bdtscore=y_train_pred[y_train == 1],weight=w_train[y_train==1])
#h.fill(dataset='Background train',bdtscore=y_train_pred[y_train == 0],weight=w_train[y_train==0])

h.scale(10)

bdtscorefig, ax = plt.subplots()

ax.set_xlabel('BDT score')
#ax.set_ylim(0.0001, 100000)
##ax.set_ylim(0.001, 10000)
ax.set_ylim(0.001, 10000)
ax.set_yscale('log')
ax.legend(title=None)
#hist.plot1d(h,overflow="over")
hist.plot1d(h)
#hist.plot1d(h,density=True)

bdtscorefig.savefig('bdtscoretest.png')

bdtscorefig, ax = plt.subplots()

ax.set_xlabel('BDT score')
#ax.set_ylim(0.0001, 100000)
##ax.set_ylim(0.001, 10000)
ax.set_ylim(0.000001, 10)
ax.set_yscale('log')
ax.legend(title=None)
#hist.plot1d(h,overflow="over")
#hist.plot1d(h)
hist.plot1d(h,density=True)

bdtscorefig.savefig('bdtscorenormtest.png')

# In[7]:


print(h.values(overflow="over"))


# In[8]:


h=coffea.hist.Hist('Events',coffea.hist.Cat('dataset', 'Dataset'),coffea.hist.Bin('bdtscore', 'BDT score', 80, 0.9, 1))
#h.fill(dataset='Signal test',bdtscore=y_test_pred[y_test == 1],weight=w_test[y_test==1])
#h.fill(dataset='Background test',bdtscore=y_test_pred[y_test == 0],weight=w_test[y_test==0])
h.fill(dataset='Signal train',bdtscore=y_train_pred[y_train == 1],weight=w_train[y_train==1])
h.fill(dataset='Background train',bdtscore=y_train_pred[y_train == 0],weight=w_train[y_train==0])

h.scale(10/9)

bdtscorefig, ax = plt.subplots()

ax.set_xlabel('BDT score')
#ax.set_ylim(0.0001, 100000)
##ax.set_ylim(0.001, 10000)
ax.set_ylim(0.001, 10000)
ax.set_yscale('log')
ax.legend(title=None)
hist.plot1d(h)
#hist.plot1d(h,density=True) 

bdtscorefig.savefig('bdtscoretrain.png')

bdtscorefig, ax = plt.subplots()

ax.set_xlabel('BDT score')
#ax.set_ylim(0.0001, 100000)
##ax.set_ylim(0.001, 10000)
ax.set_ylim(0.000001, 10)
ax.set_yscale('log')
ax.legend(title=None)
#hist.plot1d(h)
hist.plot1d(h,density=True) 

bdtscorefig.savefig('bdtscorenormtrain.png')

# In[9]:


print(h.values(overflow="over"))

