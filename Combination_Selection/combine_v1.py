# histosDict = createHisto(histosDict,"BQuarkFromLepTop",  50,  -5.,  5.)
# histosDict = createHisto(histosDict,"BQuarkFromHadTop",  50,  -5.,  5.)
# histosDict = createHisto(histosDict,"Quark0FromHadTop",  50,  -5.,  5.)
# histosDict = createHisto(histosDict,"Quark1FromHadTop",  50,  -5.,  5.)

import os
import sys
import subprocess
import argparse
import datetime
from array import array

from ROOT import *
gROOT.SetBatch(True)

sys.path.append(os.path.abspath("../../nanoAODTools"))
from treeReaderArrayTools import InputTree
from datamodel import Object

print("=========================================================")
time_start = datetime.datetime.now()
print("ProcessNanoAOD.py::START::Time("+str(time_start)+")")

path = '/root/Top_Recons/FYP2223CMS/TopReco/NanoReader/output/ntuple_topRecoTree_samples.root'
f = TFile.Open(path, 'read')
t = f.Get('topRecoTree')
#t.Show(0)
nEntriesInTree = t.GetEntries()

class ArrayObject:
    def __init__(self, pt,eta,phi,mass,JetSel_hadronFlavour,JetSel_genJetIdx, index=None):
        self._pt = pt
        self.eta=eta
        self.phi=phi
        self.mass=mass
        self.JetSel_hadronFlavour=JetSel_hadronFlavour
        self.JetSel_genJetIdx=JetSel_genJetIdx
    
def analyze(mode):
    if  isinstance(mode,int):
        SelDict = {}
        nEntriesLessFourJets = 0
        nEntriesLessFourPair = 0
        for entry in range(nEntriesInTree):
            t.GetEntry(entry)
##################################################################################################            
            if mode == 0:
                if t.nJetSel < 4:
                    nEntriesLessFourJets+=1
                    continue
                BQuarkFromLepTop = Object(t,"BQuarkFromLepTop")
                BQuarkFromHadTop = Object(t,"BQuarkFromHadTop")
                Quark0FromHadTop = Object(t,"Quark0FromHadTop")
                Quark1FromHadTop = Object(t,"Quark1FromHadTop")
                JetSel = Object(t,'JetSel')
                JetSelDict = {}
            
                for i in range(t.nJetSel):
                    JetSelDict.update({i : ArrayObject(JetSel.pt[i],JetSel.eta[i],JetSel.phi[i],JetSel.mass[i],JetSel.hadronFlavour[i],JetSel.genJetIdx[i])})
                
                entdict = {}
                for i in range(t.nJetSel):
                    if BQuarkFromLepTop.DeltaR(JetSelDict[i])<0.4:
                        entdict.update({JetSelDict[i] : BQuarkFromLepTop})
                    if BQuarkFromHadTop.DeltaR(JetSelDict[i])<0.4:
                        entdict.update({JetSelDict[i] : BQuarkFromHadTop})
                    if Quark0FromHadTop.DeltaR(JetSelDict[i])<0.4:
                        entdict.update({JetSelDict[i] : Quark0FromHadTop})
                    if Quark1FromHadTop.DeltaR(JetSelDict[i])<0.4:
                        entdict.update({JetSelDict[i] : Quark1FromHadTop})
 ##################################################################################################                       
            else:
                if t.nJetSel < mode:
                    nEntriesLessFourJets+=1
                    continue
                BQuarkFromLepTop = Object(t,"BQuarkFromLepTop")
                BQuarkFromHadTop = Object(t,"BQuarkFromHadTop")
                Quark0FromHadTop = Object(t,"Quark0FromHadTop")
                Quark1FromHadTop = Object(t,"Quark1FromHadTop")
                JetSel = Object(t,'JetSel')
                JetSelDict = {}
                
                for i in range(mode):
                    JetSelDict.update({i : ArrayObject(JetSel.pt[i],JetSel.eta[i],JetSel.phi[i],JetSel.mass[i],JetSel.hadronFlavour[i],JetSel.genJetIdx[i])})
                
                entdict = {}
                for i in range(mode):
                    if BQuarkFromLepTop.DeltaR(JetSelDict[i])<0.4:
                        entdict.update({JetSelDict[i] : BQuarkFromLepTop})
                    if BQuarkFromHadTop.DeltaR(JetSelDict[i])<0.4:
                        entdict.update({JetSelDict[i] : BQuarkFromHadTop})
                    if Quark0FromHadTop.DeltaR(JetSelDict[i])<0.4:
                        entdict.update({JetSelDict[i] : Quark0FromHadTop})
                    if Quark1FromHadTop.DeltaR(JetSelDict[i])<0.4:
                        entdict.update({JetSelDict[i] : Quark1FromHadTop})
                
            
            if len(entdict)==4:
                SelDict.update({entry : entdict})
            else : nEntriesLessFourPair +=1
        return SelDict,nEntriesLessFourJets,nEntriesLessFourPair
    else:
        print(f'Mode must be an interger or 0')
        sys.exit(1)
        


# for key, value in all.items() :
#     print (key, value)
testlist = [0,4,5,6,7,8]
for i in testlist:
    all , x, y= analyze(i)
    if i==0:
        allval=len(all)
    print(f'\n\nmode : {i}')
    print(f'Initial event : {nEntriesInTree}')
    print(f'Event with less than {i} nJetSel : {x}')
    print(f'Event with less than four paired quark {y}')
    print(f'Passed Events : {len(all)} ({len(all)/len(all)*100}%)' )
    print(f'Total events analyzed : {x+y+len(all)}')
    
time_end = datetime.datetime.now()
elapsed = time_end - time_start
elapsed_str = str(datetime.timedelta(seconds=elapsed.seconds))
print("=========================================================")
print("ProcessNanoAOD.py::DONE::Sample("+path+")::Time("+str(time_end)+")::Elapsed("+elapsed_str+")")
