# histosDict = createHisto(histosDict,"BQuarkFromLepTop",  50,  -5.,  5.)
# histosDict = createHisto(histosDict,"BQuarkFromHadTop",  50,  -5.,  5.)
# histosDict = createHisto(histosDict,"Quark0FromHadTop",  50,  -5.,  5.)
# histosDict = createHisto(histosDict,"Quark1FromHadTop",  50,  -5.,  5.)

import os
import sys
import datetime
import copy

from ROOT import *
gROOT.SetBatch(True)

sys.path.append(os.path.abspath("../../nanoAODTools"))
from datamodel import Object

print("=========================================================")
time_start = datetime.datetime.now()
print("ProcessNanoAOD.py::START::Time("+str(time_start)+")")

path = '/root/Top_Recons/FYP2223CMS/TopReco/NanoReader/outputlarge/ntuple_topRecoTree_samples.root'
file = TFile.Open(path, 'read')
tree = file.Get('topRecoTree')
#t.Show(0)
nEntriesInTree = tree.GetEntries()

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
            tree.GetEntry(entry)
            if tree.nJetSel < 4:
                nEntriesLessFourJets+=1
                continue
            BQuarkFromLepTop = Object(tree,"BQuarkFromLepTop")
            BQuarkFromHadTop = Object(tree,"BQuarkFromHadTop")
            Quark0FromHadTop = Object(tree,"Quark0FromHadTop")
            Quark1FromHadTop = Object(tree,"Quark1FromHadTop")
            JetSel = Object(tree,'JetSel')
            JetSelDict = {}
##################################################################################################            
            if tree.nJetSel<=mode or mode==0:
                for i in range(tree.nJetSel):
                    JetSelDict.update({i : ArrayObject(JetSel.pt[i],JetSel.eta[i],JetSel.phi[i],JetSel.mass[i],JetSel.hadronFlavour[i],JetSel.genJetIdx[i])})
                
                entryDict = {}
                for i in range(tree.nJetSel):
                    if BQuarkFromLepTop.DeltaR(JetSelDict[i])<0.4:
                        entryDict.update({JetSelDict[i] : BQuarkFromLepTop})
                    if BQuarkFromHadTop.DeltaR(JetSelDict[i])<0.4:
                        entryDict.update({JetSelDict[i] : BQuarkFromHadTop})
                    if Quark0FromHadTop.DeltaR(JetSelDict[i])<0.4:
                        entryDict.update({JetSelDict[i] : Quark0FromHadTop})
                    if Quark1FromHadTop.DeltaR(JetSelDict[i])<0.4:
                        entryDict.update({JetSelDict[i] : Quark1FromHadTop})
 ##################################################################################################                       
            elif tree.nJetSel>mode and mode!=0:
                for i in range(mode):
                    JetSelDict.update({i : ArrayObject(JetSel.pt[i],JetSel.eta[i],JetSel.phi[i],JetSel.mass[i],JetSel.hadronFlavour[i],JetSel.genJetIdx[i])})
                
                entryDict = {}
                for i in range(mode):
                    if BQuarkFromLepTop.DeltaR(JetSelDict[i])<0.4:
                        entryDict.update({JetSelDict[i] : BQuarkFromLepTop})
                    if BQuarkFromHadTop.DeltaR(JetSelDict[i])<0.4:
                        entryDict.update({JetSelDict[i] : BQuarkFromHadTop})
                    if Quark0FromHadTop.DeltaR(JetSelDict[i])<0.4:
                        entryDict.update({JetSelDict[i] : Quark0FromHadTop})
                    if Quark1FromHadTop.DeltaR(JetSelDict[i])<0.4:
                        entryDict.update({JetSelDict[i] : Quark1FromHadTop})
                
            
            if len(entryDict)==4:
                SelDict.update({entry : entryDict})
            else : nEntriesLessFourPair +=1
        return SelDict,nEntriesLessFourJets,nEntriesLessFourPair
    else:
        print(f'Mode must be an interger or 0')
        sys.exit(1)
        
# for key, value in all.items() :
#     print (key, value)
allval = 0
refer = 0 
testlist = [0,4,5,6,7,8,9,10,11,12,13]
with open('output.txt','w') as text:
    
    for i in testlist:
        print(f'{testlist.index(i)+1}/{len(testlist)}')
        all , x, y = analyze(i)
        if i==0: allval=len(all)
        percent = round(len(all)/allval*100,4)
        
        text.write(f'\n\nJet Combination : {i}')
        text.write(f'\nInitial event : {nEntriesInTree}')
        text.write(f'\nTotal events analyzed : {x+y+len(all)}')
        text.write(f'\nEvent with less than 4 nJetSel : {x}')
        text.write(f'\nEvent with less than four quark-jet pair {y}')
        text.write(f'\nPassed Events : {len(all)} ({percent}%)' )
        text.write(f'\nAll combination run event number : {allval}')
        if i!=0: 
            difference=round(abs(refer-percent),4)
            text.write(f'\n% diff = {difference}%')
        refer = copy.copy(percent)
    
    
time_end = datetime.datetime.now()
elapsed = time_end - time_start
elapsed_str = str(datetime.timedelta(seconds=elapsed.seconds))
print("=========================================================")
print("ProcessNanoAOD.py::DONE::Sample("+path+")::Time("+str(time_end)+")::Elapsed("+elapsed_str+")")