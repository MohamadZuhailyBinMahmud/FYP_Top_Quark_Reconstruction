import os
import sys
import datetime
import copy
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import markers
import numpy as np
from itertools import combinations,permutations
import math
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
    def __init__(self, pt,eta,phi,mass, index=None):
        self.pt = pt
        self.eta=eta
        self.phi=phi
        self.mass=mass
        
    def p4(self, corr_pt=None):
        ret = TLorentzVector()
        if corr_pt == None:
            ret.SetPtEtaPhiM(self.pt, self.eta, self.phi, self.mass)
        else:
            ret.SetPtEtaPhiM(corr_pt, self.eta, self.phi, self.mass)
        return ret

    def DeltaR(self, other):
        if isinstance(other, TLorentzVector):
            deta = abs(other.Eta() - self.eta)
            dphi = abs(other.Phi() - self.phi)
        else:
            deta = abs(other.eta - self.eta)
            dphi = abs(other.phi - self.phi)
        while dphi > math.pi:
            dphi = abs(dphi - 2 * math.pi)
        return math.sqrt(dphi**2 + deta**2)
        
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

def chisquared(perms,JetSel):
    permuationChiSave=[[],[],[],[]]
    combinationChiSve=[[],[],[],[]]
    eventChiSave=[[],[],[],[]]
    for i in perms:
        for j in i:

            #assign the jet in here
            JetBQuarkFromLepTop = ArrayObject(JetSel.pt[j[0]],JetSel.eta[j[0]],JetSel.phi[j[0]],JetSel.mass[j[0]])
            JetBQuarkFromHadTop = ArrayObject(JetSel.pt[j[1]],JetSel.eta[j[1]],JetSel.phi[j[1]],JetSel.mass[j[1]])
            JetQuark0FromHadTop = ArrayObject(JetSel.pt[j[2]],JetSel.eta[j[2]],JetSel.phi[j[2]],JetSel.mass[j[2]])
            JetQuark1FromHadTop = ArrayObject(JetSel.pt[j[3]],JetSel.eta[j[3]],JetSel.phi[j[3]],JetSel.mass[j[3]])
            
            
            #create function to calculate invariant mass
            vectorWBoson = JetQuark0FromHadTop.p4()+JetQuark1FromHadTop.p4()
            vectorHadTop = JetQuark0FromHadTop.p4()+JetQuark1FromHadTop.p4()+JetBQuarkFromHadTop.p4()
            
            #create function to calculate chisquare
            chiHadWBoson = ((80.433-vectorWBoson[3])**2)/80.433
            chiHadTop = ((172.76-vectorHadTop[3])**2)/172.76
            chiTotalHadTop = chiHadWBoson+chiHadTop
            
            #save into an array
            permuationChiSave[0].append([j[0],j[1],j[2],j[3]])
            permuationChiSave[1].append(chiHadWBoson)
            permuationChiSave[2].append(chiHadTop)
            permuationChiSave[3].append(chiTotalHadTop)
            
        
        #find the best permuation for current combination
        address = permuationChiSave[2].index(min(permuationChiSave[2]))

        #save into array
        combinationChiSve[0].append(permuationChiSave[0][address])
        combinationChiSve[1].append(permuationChiSave[1][address])
        combinationChiSve[2].append(permuationChiSave[2][address])
        combinationChiSve[3].append(permuationChiSave[3][address])
        #permuationChiSave.clear()
    
        
    #find the best combination
    address2 = combinationChiSve[2].index(min(combinationChiSve[2]))
    
    #create function to save into tree
    eventChiSave[0].append(combinationChiSve[0][address2])
    eventChiSave[1].append(combinationChiSve[1][address2])
    eventChiSave[2].append(combinationChiSve[2][address2])
    eventChiSave[3].append(combinationChiSve[3][address2])
    #combinationChiSve.clear()
    #event.update(assign(comb,eventCounter))
    #eventCounter+=1
    return eventChiSave
    
saveTree = []
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

    perms = combine(tree.nJetSel)
    savor = chisquared(perms,JetSel)
    saveTree.append(savor)

print(savor)
