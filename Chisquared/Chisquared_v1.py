import os
import sys
import datetime
import copy
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import markers
import numpy as np
from itertools import combinations,permutations

from ROOT import *
gROOT.SetBatch(True)

sys.path.append(os.path.abspath("../../../nanoAODTools"))
from datamodel import Object

print("=========================================================")
time_start = datetime.datetime.now()
print("ProcessNanoAOD.py::START::Time("+str(time_start)+")")

path = '/root/Top_Recons/FYP2223CMS/TopReco/NanoReader/output/ntuple_topRecoTree_samples.root'
file = TFile.Open(path, 'read')
tree = file.Get('topRecoTree')
nEntriesInTree = tree.GetEntries()

mode = 7

class ArrayObject:
    def __init__(self, pt,eta,phi,mass,JetSel_hadronFlavour,JetSel_genJetIdx, index=None):
        self._pt = pt
        self.eta=eta
        self.phi=phi
        self.mass=mass
        self.JetSel_hadronFlavour=JetSel_hadronFlavour
        self.JetSel_genJetIdx=JetSel_genJetIdx
        
def combine(n):
    if n<=7:
        assemble = list(range(0,n))
        
    elif n>7:
        assemble = list(range(0,7))
    
    comb = combinations(assemble,4)
    
    permarray = []
    for i in comb:
        perm = permutations(i)
        permarray.append(perm)
        
    return permarray


 
def assign(comb,eventCounter):
    jets = {}
    counter = 0
    combinatorics = {}
    
    for i in list(comb):
        combination = {'BQuarkFromLepTop':JetSelDict[i[0]],'BQuarkFromHadTop':JetSelDict[i[1]],'Quark0FromHadTop':JetSelDict[i[2]],'Quark1FromHadTop':JetSelDict[i[3]]}
        combinatorics.update({'c'+str(counter):combination})
        counter+=1
    
    jets.update({'e'+str(eventCounter):combinatorics})
    eventCounter+=1 
    return jets

event={}
for entry in range(1):
    tree.GetEntry(entry)
    if tree.nJetSel < 4:
        continue
    eventCounter = 0
    BQuarkFromLepTop = Object(tree,"BQuarkFromLepTop")
    BQuarkFromHadTop = Object(tree,"BQuarkFromHadTop")
    Quark0FromHadTop = Object(tree,"Quark0FromHadTop")
    Quark1FromHadTop = Object(tree,"Quark1FromHadTop")
    JetSel = Object(tree,'JetSel')
    JetSelDict = {}
    
    perms = combine(tree.nJetSel)
    for i in range(1):
        for j in perms[i]:
            for k in j:
                print(k)
                #assign the jet in here
                #create function to calculate invariant mass
                #create function to calculate chisquare
                #save into an array
            #find the best permuation for current combination
            #save into array
        #find the best combination
        #create function to save into tree
            
    
    #event.update(assign(comb,eventCounter))
    #eventCounter+=1
