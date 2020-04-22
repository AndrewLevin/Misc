import ROOT

for i in range(10):
    
    l1 = []

    for j in range(0,10):
        l1.append(ROOT.TH2D("","",10000,10000,0,10000,0,10000))

del l1 #if you comment this out there will be an out-of-memory crash

for i in range(10):
    
    l2 = []

    for j in range(0,10):
        l2.append(ROOT.TH2D("","",10000,10000,0,10000,0,10000))

