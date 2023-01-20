#!/usr/bin/env python3
import os
import sys
import subprocess
import argparse
import datetime
from array import array

import ROOT
ROOT.gROOT.SetBatch(True)

sys.path.append(os.path.abspath("../../nanoAODTools"))
from treeReaderArrayTools import InputTree
from datamodel import Collection
from datamodel import Event
from branchselection import BranchSelection

print("=========================================================")
time_start = datetime.datetime.now()
print("ProcessNanoAOD.py::START::Time("+str(time_start)+")")

#
# command line argument
#
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample", type=str, required=True)
args = parser.parse_args()
sample = args.sample

#
# Specify output directory
#
outDir = os.path.abspath("./output/")
# outDir = "/eos/user/n/nbinnorj/AnaFYP/TopReco/" #NOTE: Use your own eos user space

#
# Make output directory if not exist
#
if not(os.path.exists(outDir)):
  os.mkdir(outDir)
#
# Get TMPDIR
#
tmpDir = os.getenv("HOME")

#
# Read in the list of input files from txt file 
# and save them in a python list
#
inFilesList = []
txtFile = f"./small_data/{sample}.txt"
with open(txtFile) as f:
  inFilesList += f.read().splitlines()

# Remove empty lines
inFilesList = list(filter(None, inFilesList))

# Uncomment this line to test over a few number of files
# inFilesList = inFilesList[0:1] 

if len(inFilesList) == 0:
  print(f"inFilesList is empty for sample: {sample}. Please check!")
  exit()
else:
  print(f"inFilesList has {len(inFilesList)} files.")

#
# Load files into TChain
#
tChain = ROOT.TChain("Events")
for f in inFilesList:
  if f.startswith("/eos/cms/"):
    tChain.Add("root://eoscms.cern.ch/"+f)
  else:
    tChain.Add(f)

#
# Setup Tree
#
inTree  = InputTree(tChain)
nEntriesInTree = inTree.GetEntries()

#
# Enable certain branches in Events tree.
# Consider only branches that you really need for your analysis
#
branchsel = BranchSelection("branchSel.txt")
branchsel.selectBranches(inTree)

#================================
#
# Prepare histograms
#
#================================
histosDict = {}
def createHisto(hDict,name,nbin,xmin,xmax):
  hDict[name] = ROOT.TH1F(name, "", nbin, xmin, xmax)
  return hDict

histosDict = createHisto(histosDict,"h_jets_n",    15,   0, 15)
histosDict = createHisto(histosDict,"h_jets_pt",   30,   0., 300.)
histosDict = createHisto(histosDict,"h_jets_eta",  50,  -5.,  5.)
histosDict = createHisto(histosDict,"h_jets_phi",  50,  -5.,  5.)
histosDict = createHisto(histosDict,"h_jets_mass", 10,   0., 100.)
histosDict = createHisto(histosDict,"h_jets_hadronFlavour", 6, 0, 6)

#================================
#
# Place to make tree
#
#================================
makeTopRecoTree = True #values to determined whether to make root and tree

if makeTopRecoTree:
  outFileTreeName = f"ntuple_topRecoTree_{sample}.root"
  outFileTreePathTmp = tmpDir+"/"+outFileTreeName
  outFileTreePathFinal = outDir+"/"+outFileTreeName
  if outFileTreePathFinal.startswith("/eos/user"):
    outFileTreePathFinal = "root://eosuser.cern.ch/"+outFileTreePathFinal #change user name??

  outFileTree = ROOT.TFile(outFileTreePathTmp, "RECREATE")
  outFileTree.cd()
  topRecoTree = ROOT.TTree("topRecoTree", "topRecoTree")

if makeTopRecoTree: 
  def bookFloatBranch(tree, name, defaultVal): #function to allocate memory space to store branch filled with float values
    tmp = array('f', [defaultVal])
    tree.Branch(name, tmp, f'{name}/F')
    return tmp
  def bookIntBranch(tree, name, defaultVal): #function to allocate memory space to store branch filled with int values
    tmp = array('i', [defaultVal])
    tree.Branch(name, tmp, f'{name}/I')
    return tmp
  def bookFloatArrayBranch(tree, name, lenVarName, defaultN, defaultVal): #function to allocate memory space to store branch filled with float arrays
    tmp = array('f', defaultN*[defaultVal])
    tree.Branch(name, tmp, f'{name}[{lenVarName}]/F')
    return tmp
  def bookIntArrayBranch(tree, name, lenVarName, defaultN, defaultVal): #function to allocate memory space to store branch filled with int arrays
    tmp = array('i', defaultN*[defaultVal])
    tree.Branch(name, tmp, f'{name}[{lenVarName}]/I')
    return tmp

  #
  # All of these below are to create branch in the new root file
  #
  b_HadTop_pt    = bookFloatBranch(topRecoTree, 'HadTop_pt', -1.)
  b_HadTop_eta   = bookFloatBranch(topRecoTree, 'HadTop_eta', -10.)
  b_HadTop_phi   = bookFloatBranch(topRecoTree, 'HadTop_phi', -10.)
  b_HadTop_mass  = bookFloatBranch(topRecoTree, 'HadTop_mass', -1.)

  b_LepTop_pt    = bookFloatBranch(topRecoTree, 'LepTop_pt', -1.)
  b_LepTop_eta   = bookFloatBranch(topRecoTree, 'LepTop_eta', -10.)
  b_LepTop_phi   = bookFloatBranch(topRecoTree, 'LepTop_phi', -10.)
  b_LepTop_mass  = bookFloatBranch(topRecoTree, 'LepTop_mass', -1.)

  b_BQuarkFromLepTop_pt    = bookFloatBranch(topRecoTree, 'BQuarkFromLepTop_pt', -1.)
  b_BQuarkFromLepTop_eta   = bookFloatBranch(topRecoTree, 'BQuarkFromLepTop_eta', -10.)
  b_BQuarkFromLepTop_phi   = bookFloatBranch(topRecoTree, 'BQuarkFromLepTop_phi', -10.)
  b_BQuarkFromLepTop_mass  = bookFloatBranch(topRecoTree, 'BQuarkFromLepTop_mass', -1.)

  b_BQuarkFromHadTop_pt    = bookFloatBranch(topRecoTree, 'BQuarkFromHadTop_pt', -1.)
  b_BQuarkFromHadTop_eta   = bookFloatBranch(topRecoTree, 'BQuarkFromHadTop_eta', -10.)
  b_BQuarkFromHadTop_phi   = bookFloatBranch(topRecoTree, 'BQuarkFromHadTop_phi', -10.)
  b_BQuarkFromHadTop_mass  = bookFloatBranch(topRecoTree, 'BQuarkFromHadTop_mass', -1.)

  b_Quark0FromHadTop_pt    = bookFloatBranch(topRecoTree, 'Quark0FromHadTop_pt', -1.)
  b_Quark0FromHadTop_eta   = bookFloatBranch(topRecoTree, 'Quark0FromHadTop_eta', -10.)
  b_Quark0FromHadTop_phi   = bookFloatBranch(topRecoTree, 'Quark0FromHadTop_phi', -10.)
  b_Quark0FromHadTop_mass  = bookFloatBranch(topRecoTree, 'Quark0FromHadTop_mass', -1.)
  b_Quark0FromHadTop_pdgId = bookIntBranch(topRecoTree,   'Quark0FromHadTop_pdgId', -1)

  b_Quark1FromHadTop_pt    = bookFloatBranch(topRecoTree, 'Quark1FromHadTop_pt', -1.)
  b_Quark1FromHadTop_eta   = bookFloatBranch(topRecoTree, 'Quark1FromHadTop_eta', -10.)
  b_Quark1FromHadTop_phi   = bookFloatBranch(topRecoTree, 'Quark1FromHadTop_phi', -10.)
  b_Quark1FromHadTop_mass  = bookFloatBranch(topRecoTree, 'Quark1FromHadTop_mass', -1.)
  b_Quark1FromHadTop_pdgId = bookIntBranch(topRecoTree,   'Quark1FromHadTop_pdgId', -1)

  b_LepFromLepTop_pt    = bookFloatBranch(topRecoTree, 'LepFromLepTop_pt',  -1.)
  b_LepFromLepTop_eta   = bookFloatBranch(topRecoTree, 'LepFromLepTop_eta', -10.)
  b_LepFromLepTop_phi   = bookFloatBranch(topRecoTree, 'LepFromLepTop_phi', -10.)
  b_LepFromLepTop_mass  = bookFloatBranch(topRecoTree, 'LepFromLepTop_mass', -1.)
  b_LepFromLepTop_pdgId = bookIntBranch(topRecoTree,   'LepFromLepTop_pdgId', -1)

  b_NuFromLepTop_pt    = bookFloatBranch(topRecoTree, 'NuFromLepTop_pt',  -1.)
  b_NuFromLepTop_eta   = bookFloatBranch(topRecoTree, 'NuFromLepTop_eta', -10.)
  b_NuFromLepTop_phi   = bookFloatBranch(topRecoTree, 'NuFromLepTop_phi', -10.)
  b_NuFromLepTop_mass  = bookFloatBranch(topRecoTree, 'NuFromLepTop_mass', -1.)
  b_NuFromLepTop_pdgId = bookIntBranch(topRecoTree,   'NuFromLepTop_pdgId', -1)

  nMaxJets = 50
  b_nJetSel = bookIntBranch(topRecoTree, 'nJetSel', 0)
  b_JetSel_pt = bookFloatArrayBranch(topRecoTree,  'JetSel_pt',   'nJetSel', nMaxJets,  -1.)
  b_JetSel_eta = bookFloatArrayBranch(topRecoTree,  'JetSel_eta',  'nJetSel', nMaxJets,  -10.)
  b_JetSel_phi = bookFloatArrayBranch(topRecoTree,  'JetSel_phi',  'nJetSel', nMaxJets,  -10.)
  b_JetSel_mass = bookFloatArrayBranch(topRecoTree,  'JetSel_mass', 'nJetSel', nMaxJets,  -1.)
  b_JetSel_hadronFlavour = bookFloatArrayBranch(topRecoTree,'JetSel_hadronFlavour', 'nJetSel', nMaxJets, -1)
  b_JetSel_genJetIdx = bookIntArrayBranch(topRecoTree,'JetSel_genJetIdx', 'nJetSel', nMaxJets, -1)
  
  nMaxGenJets = 50
  b_nGenJet = bookIntBranch(topRecoTree, 'nGenJet', 0)
  b_GenJet_pt = bookFloatArrayBranch(topRecoTree,  'GenJet_pt',   'nGenJet', nMaxGenJets,  -1.)
  b_GenJet_eta = bookFloatArrayBranch(topRecoTree,  'GenJet_eta',  'nGenJet', nMaxGenJets,  -10.)
  b_GenJet_phi = bookFloatArrayBranch(topRecoTree,  'GenJet_phi',  'nGenJet', nMaxGenJets,  -10.)
  b_GenJet_mass = bookFloatArrayBranch(topRecoTree,  'GenJet_mass', 'nGenJet', nMaxGenJets,  -1.)
  b_GenJet_hadronFlavour = bookIntArrayBranch(topRecoTree,'GenJet_hadronFlavour', 'nGenJet', nMaxGenJets, -1)

#
# https://github.com/scikit-hep/particle/blob/master/src/particle/data/particle2022.csv
#
dict_pdgId_particlemass = {
  1: 0.00467,#up
  2: 0.00216,#down
  3: 0.0934,#strange
  4: 1.270,#charm
  5: 4.180,#bottom
  11: 0.00051099895,#electron
  12: 0.,#electron neutrino
  13: 0.1056583755,#muon 
  14: 0.,#muon neutrino
  15: 1.77686,#tau
  16: 0.,#tau neutrino
}
#================================
#
# Start the event loop
#
#================================
nMaxEntries = 100000

nEntriesToLoop = nEntriesInTree
if (nMaxEntries > 0) and nEntriesInTree > nMaxEntries:
  nEntriesToLoop = nMaxEntries 

print(f"nEntriesInTree={nEntriesInTree}")
print(f"nEntriesToLoop={nEntriesToLoop}")

for iEntry in range(0, nEntriesToLoop):
  if iEntry%1000 == 0:
    print(f"{iEntry} / {nEntriesToLoop}") #??

  ###########################################
  #
  # Load the event from the tree
  #
  ###########################################
  event = Event(inTree, iEntry)

  ###########################################
  #
  # Load the event from the tree
  #
  ###########################################
  b_nGenJet[0] = 0 #??
  b_nJetSel[0] = 0 #??

  ############################################
  #
  # Get all jets from jet collection in tree
  #
  ###########################################
  def loadGenHistory(event):
    try:
      genparts = event.genparts
    except RuntimeError as e:
      genparts = Collection(event, "GenPart") #this is the supposed oop-like module
      for idx, gp in enumerate(genparts):
        if 'dauIdx' not in gp.__dict__:
          gp.dauIdx = []
        if gp.genPartIdxMother >= 0:
          mom = genparts[gp.genPartIdxMother]
          if 'dauIdx' not in mom.__dict__:
            mom.dauIdx = [idx]
          else:
            mom.dauIdx.append(idx)
      event.genparts = genparts
    #
    #
    #
    def isHadronic(gp):
      if len(gp.dauIdx) == 0:
        raise ValueError('isHadronic():Particle has no daughters!')
      for idx in gp.dauIdx:
        if abs(genparts[idx].pdgId) < 6:
          return True
      return False
    def isLeptonic(gp):
      if len(gp.dauIdx) == 0:
        raise ValueError('isLeptonic():Particle has no daughters!')
      for idx in gp.dauIdx:
        if abs(genparts[idx].pdgId) >= 11 and abs(genparts[idx].pdgId) <= 16:
          return True
      return False
    def getFinal(gp):
      for idx in gp.dauIdx:
        dau = genparts[idx]
        if dau.pdgId == gp.pdgId:
          return getFinal(dau)
        return gp

    event.nGenTops = 0
    event.nGenWs = 0

    event.GenTopsForMtt = []
    event.GenAntiTopsForMtt = []

    event.hadGenTops = []
    event.lepGenTops = []

    event.hadGenWs = []
    event.lepGenWs = []

    for gp in genparts:
      #
      # Only for Pow+Pythia
      #
      if (gp.pdgId == 6) and (gp.status==22):
        event.GenTopsForMtt.append(gp)
      elif (gp.pdgId == -6) and (gp.status==22):
        event.GenAntiTopsForMtt.append(gp)
      #
      #
      #
      if gp.statusFlags & (1 << 13) == 0:
        continue
      #
      #
      #
      if abs(gp.pdgId) == 6:
        event.nGenTops += 1
        for idx in gp.dauIdx:
          dau = genparts[idx]
          if abs(dau.pdgId) == 24:
            genW = getFinal(dau)
            gp.genW = genW
            if isHadronic(genW):
              event.hadGenTops.append(gp)
            if isLeptonic(genW):
              event.lepGenTops.append(gp)
          elif abs(dau.pdgId) in (1, 3, 5):
            gp.genBQuark = dau
      #
      # 
      #
      elif abs(gp.pdgId) == 24:
        event.nGenWs += 1
        if isHadronic(gp):
          event.hadGenWs.append(gp)
        if isLeptonic(gp):
          event.lepGenWs.append(gp)
    event.genparts = genparts

  loadGenHistory(event)
  #
  # get the tops
  #
  lepTop = event.lepGenTops[0]
  hadTop = event.hadGenTops[0]

  #
  # get the quarks
  #
  bquark_lepTop = lepTop.genBQuark
  bquark_hadTop = hadTop.genBQuark
  hadWFromTop = event.hadGenWs[0]
  # Sort by pt
  if event.genparts[hadWFromTop.dauIdx[0]].pt > event.genparts[hadWFromTop.dauIdx[1]].pt:
    quark0_hadWFromTop = event.genparts[hadWFromTop.dauIdx[0]]
    quark1_hadWFromTop = event.genparts[hadWFromTop.dauIdx[1]]
  else:
    quark0_hadWFromTop = event.genparts[hadWFromTop.dauIdx[1]]
    quark1_hadWFromTop = event.genparts[hadWFromTop.dauIdx[0]]

  #
  # get the lepton and neutrino
  #
  lepWFromTop = event.lepGenWs[0]

  if abs(event.genparts[lepWFromTop.dauIdx[0]].pdgId) in [11, 13, 15]:
    lep_lepWFromTop = event.genparts[lepWFromTop.dauIdx[0]]
    nu_lepWFromTop  = event.genparts[lepWFromTop.dauIdx[1]]
  else:
    lep_lepWFromTop = event.genparts[lepWFromTop.dauIdx[1]]
    nu_lepWFromTop  = event.genparts[lepWFromTop.dauIdx[0]]
  
  #
  # Lets only consider real electron or muon channel ttbar events
  #
  if not(abs(lep_lepWFromTop.pdgId)==11 or abs(lep_lepWFromTop.pdgId)==13): continue

  b_HadTop_pt[0]   = hadTop.pt
  b_HadTop_eta[0]  = hadTop.eta
  b_HadTop_phi[0]  = hadTop.phi
  b_HadTop_mass[0] = hadTop.mass

  b_LepTop_pt[0]   = lepTop.pt
  b_LepTop_eta[0]  = lepTop.eta
  b_LepTop_phi[0]  = lepTop.phi
  b_LepTop_mass[0] = lepTop.mass

  b_BQuarkFromLepTop_pt[0]    = bquark_lepTop.pt
  b_BQuarkFromLepTop_eta[0]   = bquark_lepTop.eta
  b_BQuarkFromLepTop_phi[0]   = bquark_lepTop.phi
  b_BQuarkFromLepTop_mass[0]  = dict_pdgId_particlemass[abs(bquark_lepTop.pdgId)]

  b_BQuarkFromHadTop_pt[0]    = bquark_hadTop.pt
  b_BQuarkFromHadTop_eta[0]   = bquark_hadTop.eta
  b_BQuarkFromHadTop_phi[0]   = bquark_hadTop.phi
  b_BQuarkFromHadTop_mass[0]  = dict_pdgId_particlemass[abs(bquark_hadTop.pdgId)]

  b_Quark0FromHadTop_pt[0]    = quark0_hadWFromTop.pt
  b_Quark0FromHadTop_eta[0]   = quark0_hadWFromTop.eta
  b_Quark0FromHadTop_phi[0]   = quark0_hadWFromTop.phi
  b_Quark0FromHadTop_mass[0]  = dict_pdgId_particlemass[abs(quark0_hadWFromTop.pdgId)]
  b_Quark0FromHadTop_pdgId[0] = quark0_hadWFromTop.pdgId

  b_Quark1FromHadTop_pt[0]    = quark1_hadWFromTop.pt
  b_Quark1FromHadTop_eta[0]   = quark1_hadWFromTop.eta
  b_Quark1FromHadTop_phi[0]   = quark1_hadWFromTop.phi
  b_Quark1FromHadTop_mass[0]  = dict_pdgId_particlemass[abs(quark1_hadWFromTop.pdgId)]
  b_Quark1FromHadTop_pdgId[0] = quark1_hadWFromTop.pdgId

  b_LepFromLepTop_pt[0]    = lep_lepWFromTop.pt
  b_LepFromLepTop_eta[0]   = lep_lepWFromTop.eta
  b_LepFromLepTop_phi[0]   = lep_lepWFromTop.phi
  b_LepFromLepTop_mass[0]  = dict_pdgId_particlemass[abs(lep_lepWFromTop.pdgId)]
  b_LepFromLepTop_pdgId[0] = lep_lepWFromTop.pdgId

  b_NuFromLepTop_pt[0]    = nu_lepWFromTop.pt
  b_NuFromLepTop_eta[0]   = nu_lepWFromTop.eta
  b_NuFromLepTop_phi[0]   = nu_lepWFromTop.phi
  b_NuFromLepTop_mass[0]  = nu_lepWFromTop.mass
  b_NuFromLepTop_pdgId[0] = nu_lepWFromTop.pdgId

  ############################################
  #
  # Get all jets from jet collection in tree
  #
  ###########################################
  genJetsAll = Collection(event, "GenJet")

  b_nGenJet[0] = len(genJetsAll)
  for iJet, jet in enumerate(genJetsAll): 
    if makeTopRecoTree:
      b_GenJet_pt[iJet] = jet.pt
      b_GenJet_eta[iJet] = jet.eta
      b_GenJet_phi[iJet] = jet.phi
      b_GenJet_mass[iJet] = jet.mass
      b_GenJet_hadronFlavour[iJet] = jet.hadronFlavour

  jetsAll = Collection(event, "Jet")
  #
  # Select jets
  #
  jetsSel = [x for x in jetsAll 
    if abs(x.eta) < 2.4 
    and x.pt > 20.
    and x.jetId & (1<<1) # Tight JetID
    and x.puId & (1<<2)# Pileup Jet ID Loose
  ]
  
  b_nJetSel[0] = len(jetsSel)
  for iJet, jet in enumerate(jetsSel): 
    #
    #
    #
    histosDict["h_jets_pt"].Fill(jet.pt)
    histosDict["h_jets_eta"].Fill(jet.eta)
    histosDict["h_jets_phi"].Fill(jet.phi)
    histosDict["h_jets_mass"].Fill(jet.mass)
    histosDict["h_jets_hadronFlavour"].Fill(jet.hadronFlavour)
    #
    #
    #
    if makeTopRecoTree:
      b_JetSel_pt[iJet] = jet.pt
      b_JetSel_eta[iJet] = jet.eta
      b_JetSel_phi[iJet] = jet.phi
      b_JetSel_mass[iJet] = jet.mass
      b_JetSel_hadronFlavour[iJet] = jet.hadronFlavour
      b_JetSel_genJetIdx[iJet] = jet.genJetIdx

  #
  #
  #
  if makeTopRecoTree:
    topRecoTree.Fill()
#================================
#
# End the event loop
#
#================================

#
# Close outFileTree 
#
if makeTopRecoTree:
  outFileTree.cd()
  topRecoTree.Write()
  outFileTree.Close()

  print("=========================================================")
  print("Copying tree file to final destination directory: "+outDir)
  if outFileTreePathFinal.startswith("root://"):
    command = ["xrdcp", "-f", outFileTreePathTmp, outFileTreePathFinal]
  else:
    command = ["rsync", "-v", outFileTreePathTmp, outFileTreePathFinal]
  subprocess.call(command)
  print("=========================================================")
  print("Deleting local temporary tree file: "+outFileTreePathTmp)
  command = ["rm", "-fv", outFileTreePathTmp]
  subprocess.call(command) #call() is a py2 command

#
# Save all histograms in outFileHisto
#
outFileHistoName = f"Histos_{sample}.root"
outFileHistoPathFinal = outDir+"/"+outFileHistoName
if outFileHistoPathFinal.startswith("/eos/user"):
  outFileHistoPathFinal = "root://eosuser.cern.ch/"+outFileHistoPathFinal

outFileHisto = ROOT.TFile(outFileHistoPathFinal,"RECREATE")
outFileHisto.cd()
for hName in histosDict:
  histosDict[hName].Write()
outFileHisto.Close()

time_end = datetime.datetime.now()
elapsed = time_end - time_start
elapsed_str = str(datetime.timedelta(seconds=elapsed.seconds))
print("=========================================================")
print("ProcessNanoAOD.py::DONE::Sample("+sample+")::Time("+str(time_end)+")::Elapsed("+elapsed_str+")")
