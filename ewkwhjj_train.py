import torch
from coffea.util import load

from coffea import hist

import coffea.hist

import pandas

from sklearn.metrics import roc_curve, auc

import matplotlib.pyplot as plt

import mplhep as hep

result_2018 = load("training_data_2018")

resolved_variables = ['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','met','metphi','higgsjet1pt','higgsjet2pt','vbfjet1pt','vbfjet2pt','higgsjet1eta','higgsjet2eta','vbfjet1eta','vbfjet2eta','higgsjet1phi','higgsjet2phi','vbfjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbfjet1btag','vbjfet2btag','higgsdijetmass','vbfdijetmass','vbfdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt','higgdijetdeltar','vbfdijetdeltar','higgsdijetabsdeta','dr1','dr2','dr3','dr4']


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

training_data = training_data_resolved_2018

result = result_2018

training_data = training_data.reset_index(drop=True)

lumi=59.6

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

for k in xs.keys():
    if 'wlepht' in k:
        xs[k] = xs[k]

training_data = training_data[(training_data["label"] == "ewkwhjj")
                              | (training_data["label"] == "ttsemi")
                              | (training_data["label"] == "stoptchan")
                              | (training_data["label"] == "santitoptchan")
                              | training_data["label"].str.contains("wlepht")
                                 ]

training_data= training_data.reset_index(drop=True)

labels = []
weights = []

for i in range(len(training_data)):
    assert(training_data["label"][i] in xs.keys())    

    weights.append(xs[training_data["label"][i]]*1000*lumi/result["nevents"][training_data["label"][i]])
        
    if training_data["label"][i] == "ewkwhjj":
        labels.append(0)
    elif training_data["label"][i] == "ttsemi":
        labels.append(1)
    elif training_data["label"][i] == "stoptchan":
        labels.append(2)
    elif training_data["label"][i] == "santitoptchan":
        labels.append(3)        
    else:    
        labels.append(4)
        
assert(len(weights) == len(labels))
assert(len(weights) == len(training_data))

labels=pandas.DataFrame(labels,columns=["label"])
weights=pandas.DataFrame(weights,columns=["weight"])

training_data = training_data[training_data.columns.drop(["weight","label"])]

training_data = pandas.concat([training_data,labels,weights],axis=1)

print(len(training_data[training_data["label"] == 0]))

print(len(training_data[training_data["label"] == 1]))

X = training_data[training_data.columns.drop(["label","weight"])]

w = training_data["weight"]

from sklearn.model_selection import train_test_split 

y = training_data["label"]

X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(X, y, w, test_size=0.5, random_state=1)

print(len(X_train))
print(len(X_test))

#print(X_train)
#print(X_test)

#print(X_train.shape)

import GPUtil
print(GPUtil.getAvailable())

use_cuda = torch.cuda.is_available()

#use_cuda = False

device = torch.device("cuda:2" if use_cuda else "cpu")

print("Device: ",device)

print('__CUDNN VERSION:', torch.backends.cudnn.version())
print('__Number CUDA Devices:', torch.cuda.device_count())
print('__CUDA Device Name:',torch.cuda.get_device_name(1))
print('__CUDA Device Total Memory [GB]:',torch.cuda.get_device_properties(1).total_memory/1e9)

from torch import nn
from torch.utils.data import DataLoader
from torchvision import datasets, transforms

class NeuralNetwork(nn.Module):
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(37, 512),
            nn.ReLU(),
#            nn.Sigmoid(),
            nn.Linear(512, 512),
            nn.ReLU(),
#            nn.Linear(512, 512),
#            nn.ReLU(),
#            nn.Linear(512, 512),
#            nn.ReLU(),
#            nn.Linear(512, 512),
#            nn.ReLU(),
#            nn.Linear(512, 512),
#            nn.ReLU(),            
#            nn.Sigmoid(),
#            nn.Linear(512, 512),
#            nn.ReLU(),                        
#            nn.Sigmoid(),
            nn.Linear(512, 5),
#            nn.ReLU(),                        
        )

    def forward(self, x):
#        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits

model = NeuralNetwork().to(device)
print(model)

#loss_fn = nn.CrossEntropyLoss()
#loss_fn = nn.CrossEntropyLoss(weight=torch.tensor([1.,1.,1.,1.,1.],device=device))
bg_tot_xs = xs['ttsemi']+xs['stoptchan']+xs['santitoptchan']+2000.
#loss_fn = nn.CrossEntropyLoss(weight=torch.tensor([1.,xs['ttsemi']/bg_tot_xs,xs['stoptchan']/bg_tot_xs,xs['santitoptchan']/bg_tot_xs,2000./bg_tot_xs],device=device))
loss_fn = nn.CrossEntropyLoss(weight=torch.tensor([1.,1.,1.,1.,1.],device=device))


#loss_fn = nn.NLLLoss()
#loss_fn = nn.MSELoss() 

learning_rate = 1e-3
#batch_size = 64
batch_size = 256

#optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

import numpy as np



#t1 = torch.from_numpy(X_train.to_numpy(),dtype=torch.float)
t1 = torch.tensor(torch.from_numpy(np.array(X_train.to_numpy(),dtype="f")),device=device)


#print(list(torch.utils.data.RandomSampler(t1, replacement=False)))
#print(t1[list(torch.utils.data.RandomSampler(t1, replacement=False))])

t1_test = torch.tensor(torch.from_numpy(np.array(X_test.to_numpy(),dtype="f")),device=device)

#t2 = torch.tensor(torch.from_numpy(np.array(np.transpose([y_train.to_numpy() == 0,y_train.to_numpy() == 1, y_test.to_numpy() == 2]),dtype="f")),device=device)
#t2_test = torch.tensor(torch.from_numpy(np.array(np.transpose([y_test.to_numpy() == 0,y_test.to_numpy() == 1, y_test.to_numpy() == 2]),dtype="f")),device=device)

#t2 = torch.tensor(torch.from_numpy(np.array(np.transpose([y_train.to_numpy(),1-y_train.to_numpy()]),dtype="f")),device=device)
#t2_test = torch.tensor(torch.from_numpy(np.array(np.transpose([y_test.to_numpy(),1-y_test.to_numpy()]),dtype="f")),device=device)

#t2 = torch.tensor(torch.from_numpy(np.array(np.transpose([y_train.to_numpy(),y_train.to_numpy()]),dtype="f")),device=device)
#t2_test = torch.tensor(torch.from_numpy(np.array(np.transpose([y_test.to_numpy(),y_test.to_numpy()]),dtype="f")),device=device)

t2 = torch.tensor(torch.from_numpy(np.array(y_train.to_numpy(),dtype="uint8")),device=device)
t2_test = torch.tensor(torch.from_numpy(np.array(y_test.to_numpy(),dtype="uint8")),device=device)

#t2 = torch.tensor(torch.from_numpy(np.array(np.transpose([y_train.to_numpy(),1-y_train.to_numpy()]),dtype="l")),device=device)
#t2 = torch.from_numpy(np.array(y_train.to_numpy(),dtype="l"))

print(t1.dtype)
print(t2.dtype)
print(t1_test.dtype)
print(t2_test.dtype)
print(t1.size())
print(t2.size())
print(t1_test.size())
print(t2_test.size())



epoch_list = []
train_loss_list = []
test_loss_list = []

#print(list(torch.utils.data.RandomSampler(t1, replacement=False)))
#print(list(torch.utils.data.BatchSampler(torch.utils.data.RandomSampler(t1, replacement=False),batch_size,True)))

#epochs = 50
#epochs = 1000
epochs = 10
for e in range(epochs):
    if e % 1 == 0:
        print(e)

    epoch_list.append(e)    

    for batch_indices in list(torch.utils.data.BatchSampler(torch.utils.data.RandomSampler(t1, replacement=False),batch_size,True)):

        optimizer.zero_grad()

        #shuffled_indices = list(torch.utils.data.RandomSampler(t1, replacement=False))
    
        t1_batch = t1[batch_indices]
        t2_batch = t2[batch_indices]

        output = model(t1_batch)

        loss = loss_fn(output, t2_batch)

        loss.backward()        

        optimizer.step()

    with torch.no_grad():
        
        output_train = model(t1)
        
        train_loss_list.append(loss_fn(output_train, t2).item())

        output_test =  nn.Softmax(dim=1)(model(t1_test))

        test_loss_list.append(loss_fn(output_test, t2_test).item())
    
        if e % 1 == 0:
            print(loss.item())
        
            output_train =  nn.Softmax(dim=1)(model(t1))

            print(loss_fn(output_test, t2_test).item())
        
            output_train = output_train.cpu()
            output_train = output_train.detach()
            y_train_pred = output_train.numpy()[:,0]
            train_false_positive_rates, train_true_positive_rates, _ = roc_curve(np.array(y_train == 0,dtype=int), y_train_pred[:])
            print(auc(train_false_positive_rates, train_true_positive_rates))
    
            output_test = output_test.cpu()
            output_test = output_test.detach()
            y_test_pred = output_test.numpy()[:,0]
            test_false_positive_rates, test_true_positive_rates, _ = roc_curve(np.array(y_test == 0,dtype=int), y_test_pred[:])
            print(auc(test_false_positive_rates, test_true_positive_rates))

torch.save(model.state_dict(), 'model_weights.pth')

output_train = output_train.cpu()
output_train = output_train.detach()
y_train_pred = output_train.numpy()[:,0]
#train_false_positive_rates, train_true_positive_rates, _ = roc_curve(y_train[:], y_train_pred[:])
train_false_positive_rates, train_true_positive_rates, _ = roc_curve(np.array(y_train == 0,dtype=int), y_train_pred[:])


print(output_test)
output_test_softmax = nn.Softmax(dim=1)(output_test)
print(output_test_softmax)

print(auc(train_false_positive_rates, train_true_positive_rates))

output_test = output_test.cpu()
output_test = output_test.detach()
y_test_pred = output_test.numpy()[:,0]
test_false_positive_rates, test_true_positive_rates, _ = roc_curve(np.array(y_test == 0,dtype=int), y_test_pred[:])
print(auc(test_false_positive_rates, test_true_positive_rates))


output_test_softmax = output_test_softmax.cpu()
output_test_softmax = output_test_softmax.detach()
y_test_softmax_pred = output_test_softmax.numpy()[:,0]
test_softmax_false_positive_rates, test_softmax_true_positive_rates, _ = roc_curve(np.array(y_test == 0,dtype=int), y_test_softmax_pred[:])
print(auc(test_softmax_false_positive_rates, test_softmax_true_positive_rates))

plt.style.use(hep.style.ROOT) 

roc_fig = plt.figure()

plt.plot(1-test_false_positive_rates, test_true_positive_rates)
plt.plot(1-train_false_positive_rates, train_true_positive_rates)

plt.xlabel('Background Rejection Rate')
plt.ylabel('Signal Efficiency')

roc_fig.savefig("roc.png")

loss_fig = plt.figure()

plt.plot(np.array(epoch_list), np.array(train_loss_list))
plt.plot(np.array(epoch_list), np.array(test_loss_list))

plt.xlabel('Epoch')
plt.ylabel('Cross Entropy Loss')

loss_fig.savefig("loss.png")

#h=coffea.hist.Hist('Events',coffea.hist.Cat('dataset', 'Dataset'),coffea.hist.Bin('bdtscore', 'BDT score', 100, 0, 1))
h=coffea.hist.Hist('Events',coffea.hist.Cat('dataset', 'Dataset'),coffea.hist.Bin('bdtscore', 'BDT score', 50, 0, 1))

print(y_train_pred)

h.fill(dataset='Signal train',bdtscore=y_train_pred[y_train == 1],weight=w_train[y_train==1])
h.fill(dataset='Background train',bdtscore=y_train_pred[y_train == 0],weight=w_train[y_train==0])
#h.fill(dataset='Signal train',bdtscore=y_train_pred[y_train == 1])
#h.fill(dataset='Background train',bdtscore=y_train_pred[y_train == 0])

print(np.sum(h.values(overflow='all')[('Background train',)]))
print(result["nevents"]['wlepsherpa'])

bdtscorefig, ax = plt.subplots()

ax.set_xlabel('BDT score')
ax.set_ylim(0.001, 10000)
ax.set_yscale('log')
ax.legend(title=None)

hist.plot1d(h,overflow="over")

bdtscorefig.savefig('bdtscore.png')
