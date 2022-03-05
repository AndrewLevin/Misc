import torch

from coffea.util import load

import pandas

from sklearn.metrics import roc_curve, auc

import matplotlib.pyplot as plt

import mplhep as hep

training_data = load("ewkwhjj_training_data")

print(len(training_data[training_data["label"] == 0]))

print(len(training_data[training_data["label"] == 1]))

X = training_data[training_data.columns.drop(["label","weight"])]

from sklearn.model_selection import train_test_split 

y = training_data["label"]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=1)

print(len(X_train))
print(len(X_test))

#print(X_train)
#print(X_test)

#print(X_train.shape)

import GPUtil
print(GPUtil.getAvailable())

use_cuda = torch.cuda.is_available()

#use_cuda = False

device = torch.device("cuda" if use_cuda else "cpu")

print("Device: ",device)

print('__CUDNN VERSION:', torch.backends.cudnn.version())
print('__Number CUDA Devices:', torch.cuda.device_count())
print('__CUDA Device Name:',torch.cuda.get_device_name(0))
print('__CUDA Device Total Memory [GB]:',torch.cuda.get_device_properties(0).total_memory/1e9)

from torch import nn
from torch.utils.data import DataLoader
from torchvision import datasets, transforms

class NeuralNetwork(nn.Module):
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(30, 512),
            nn.ReLU(),
#            nn.Sigmoid(),
            nn.Linear(512, 512),
            nn.ReLU(),
            nn.Linear(512, 512),
            nn.ReLU(),
            nn.Linear(512, 512),
            nn.ReLU(),
            nn.Linear(512, 512),
            nn.ReLU(),
            nn.Linear(512, 512),
            nn.ReLU(),            
#            nn.Sigmoid(),
            nn.Linear(512, 512),
            nn.ReLU(),                        
#            nn.Sigmoid(),
            nn.Linear(512, 2),
#            nn.ReLU(),                        
#            nn.Sigmoid()
        )

    def forward(self, x):
#        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits

model = NeuralNetwork().to(device)
print(model)

loss_fn = nn.CrossEntropyLoss()
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

t2 = torch.tensor(torch.from_numpy(np.array(np.transpose([y_train.to_numpy(),1-y_train.to_numpy()]),dtype="f")),device=device)
t2_test = torch.tensor(torch.from_numpy(np.array(np.transpose([y_test.to_numpy(),1-y_test.to_numpy()]),dtype="f")),device=device)
#t2 = torch.tensor(torch.from_numpy(np.array(np.transpose([y_train.to_numpy(),1-y_train.to_numpy()]),dtype="l")),device=device)
#t2 = torch.from_numpy(np.array(y_train.to_numpy(),dtype="l"))

epoch_list = []
train_loss_list = []
test_loss_list = []

#print(list(torch.utils.data.RandomSampler(t1, replacement=False)))
#print(list(torch.utils.data.BatchSampler(torch.utils.data.RandomSampler(t1, replacement=False),batch_size,True)))

#epochs = 50
epochs = 1000
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
#    output = model(t1)
#    output = nn.Softmax(dim=1)(model(t1))

#    loss = loss_fn(output[:,0], t2[:,0])
        loss = loss_fn(output, t2_batch)

        loss.backward()        

        optimizer.step()

    output_train = model(t1)
        
    train_loss_list.append(loss_fn(output_train, t2).item())

    output_test =  model(t1_test)

    test_loss_list.append(loss_fn(output_test, t2_test).item())
    
    if e % 1 == 0:
#        print(output)
#        print(t2)
        print(loss.item())
        
#        output_test =  nn.Softmax(dim=1)(model(t1_test))

        output_train =  model(t1)

        print(loss_fn(output_test, t2_test).item())
        
        output_train = output_train.cpu()
        output_train = output_train.detach()
        y_train_pred = output_train.numpy()[:,0]
        train_false_positive_rates, train_true_positive_rates, _ = roc_curve(y_train[:], y_train_pred[:])
        print(auc(train_false_positive_rates, train_true_positive_rates))
    
        output_test = output_test.cpu()
        output_test = output_test.detach()
        y_test_pred = output_test.numpy()[:,0]
        test_false_positive_rates, test_true_positive_rates, _ = roc_curve(y_test[:], y_test_pred[:])
        print(auc(test_false_positive_rates, test_true_positive_rates))

torch.save(model.state_dict(), 'model_weights.pth')

output_train = output_train.cpu()
output_train = output_train.detach()
y_train_pred = output_train.numpy()[:,0]
train_false_positive_rates, train_true_positive_rates, _ = roc_curve(y_train[:], y_train_pred[:])


print(output_test)
output_test_softmax = nn.Softmax(dim=1)(output_test)
print(output_test_softmax)

print(auc(train_false_positive_rates, train_true_positive_rates))

output_test = output_test.cpu()
output_test = output_test.detach()
y_test_pred = output_test.numpy()[:,0]
test_false_positive_rates, test_true_positive_rates, _ = roc_curve(y_test[:], y_test_pred[:])
print(auc(test_false_positive_rates, test_true_positive_rates))


output_test_softmax = output_test_softmax.cpu()
output_test_softmax = output_test_softmax.detach()
y_test_softmax_pred = output_test_softmax.numpy()[:,0]
test_softmax_false_positive_rates, test_softmax_true_positive_rates, _ = roc_curve(y_test[:], y_test_softmax_pred[:])
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
