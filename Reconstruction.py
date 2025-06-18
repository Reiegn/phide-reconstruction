import ROOT
from itertools import permutations
from concurrent.futures import ProcessPoolExecutor
from time import time
import math

start_time = time()

ROOT.gSystem.Load("libDelphes")

try:
  ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
  ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
  pass


# Global chi_squared function
def chi_squared(candidate):
    W1 = (candidate[0].P4() + candidate[1].P4()).M()
    W2 = (candidate[3].P4() + candidate[4].P4()).M()
    t1 = (candidate[0].P4() + candidate[1].P4() + candidate[2].P4()).M()
    t2 = (candidate[3].P4() + candidate[4].P4() + candidate[5].P4()).M()
    phi = (candidate[6].P4() + candidate[7].P4()).M()
    
    #Standard Deviation of same reconstruction with unweighted chisq.
    sw, st, sp = 17.07, 36.91, 21.46 #21.46, 42.92, 64.38, 85.84, 106.92, 128.76

    chisq = ((W1 - 81.3) ** 2 / sw ** 2) + ((W2 - 81.3) ** 2 / sw ** 2) + ((t1 - 172.6) ** 2 / st ** 2) + ((t2 - 172.6) ** 2 / st ** 2) #+ ((phi - 100) ** 2 / sp ** 2)

    return chisq, W1, W2, t1, t2, phi

# Process a single chunk to find the best permutation.
def process_chunk(chunk, chi_squared):
    best_permutation = min(chunk, key=lambda x: chi_squared(x)[0])
    return best_permutation, chi_squared(best_permutation)

#Divides an events possible permutations into chunks.
def chunkify(permutations, num_chunks):

    chunk_size = math.ceil(len(permutations) / num_chunks)

    return [permutations[i:i + chunk_size] for i in range(0, len(permutations), chunk_size)]

# 1 hadronic decay and 1 leptonic decay reconstruction
def half_reconstruction(selectedJets, lepton, MET, max_workers):

    output = []
    output = [(j[0], j[1], j[2], lepton, MET, j[3], j[4], j[5])
              for j in permutations(selectedJets, 6)]

    #'output' divided into chunks
    chunks = chunkify(output, max_workers)

    #parallel computation of chunks
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        
        #Runs process_chunk for each chunk, creating a list of each chunk's permutation with lowest chisq.
        results = executor.map(process_chunk, chunks, [chi_squared] * len(chunks))

    best_permutation, best_values = min(results, key=lambda res: res[1][0])
    _, W1, W2, t1, t2, phi = best_values

    return W1, W2, t1, t2, phi

# full_hadronic reconstruction
def full_hadronic_reconstruction(selectedJets, max_workers):

    output = []
    output = [(j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7])
              for j in permutations(selectedJets, 8)]

    chunks = chunkify(output, max_workers)

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(process_chunk, chunks, [chi_squared] * len(chunks))

    best_permutation, best_values = min(results, key=lambda res: res[1][0])
    _, W1, W2, t1, t2, phi = best_values

    return W1, W2, t1, t2, phi

# Returns filtered list by highest PT, cutoff at a length of 2 more than the minimum required for event.
def sort_jets_by_pt(jets, count):
    sorted_jets = sorted(jets, key=lambda jet: jet.PT, reverse=True)
    
    best_jets = sorted_jets[:count+1]
    return best_jets

# Runs reconstructions based off event's decay mode.
def execute(category, HW, Ht, Hphi, selectedJets, lepton, MET, max_workers=4):
    if category == 1:
        best_jets = sort_jets_by_pt(selectedJets, 8)
        W1, W2, t1, t2, phi = full_hadronic_reconstruction(best_jets, max_workers)
        HW.Fill(W1)
        HW.Fill(W2)
        Ht.Fill(t1)
        Ht.Fill(t2)
        Hphi.Fill(phi)

    elif category == 3:
        best_jets = sort_jets_by_pt(selectedJets, 6)
        W1, W2, t1, t2, phi = half_reconstruction(best_jets, lepton[0], MET, max_workers)
        HW.Fill(W1)
        HW.Fill(W2)
        Ht.Fill(t1)
        Ht.Fill(t2)
        Hphi.Fill(phi)

# Strict event selection for all possibilities of dual-top decay, accounting for extra 2 jets from phide.
def decay_channels(l):
  
 # l = [len(leptons), len(selectedJets), len(ljets), len(bjets)]

  #Not enough jets detected for any decay mode.
  if l[1] < 4: return 0

  #0 Leptons - Hadronic
  elif l[0] == 0 and l[1] >= 8 and l[2] >= 4: return 1

  #1 Lepton
  elif l[0] == 1:
    
    #Hadronic or Half
    if l[1] >= 8 and l[2] >= 4: return 2
    
    #Half
    elif l[1] >= 6 and l[2] >= 2: return 3

    else: return 7

  #At least 2 Leptons
  elif l[0] >= 2:

    #All
    if l[1] >= 8 and l[2] >= 4: return 4

    #Half or Leptonic
    elif l[1] >= 6 and l[2] >= 2: return 5

    #Leptonic
    elif l[1] >= 4: return 6

    else: return 7

  #Not enough detected products for reconstruction.
  else: return 7

def ljetfilter(jets):
    def jetfilter(jet): return jet.BTag==0
    return filter(jetfilter,jets)

def bjetfilter(jets):
  def jetfilter(jet): return jet.BTag==1
  return filter(jetfilter,jets)

def jetSelection(jets, ptcut, etacut):
    def jetfilter(jet): return jet.PT>ptcut and abs(jet.Eta)<etacut
    return filter(jetfilter,jets)
    
def analyze(start, end, tree, branchJet, branchElectron, branchMuon, branchMET, HW, Ht, Hphi):
    
    for entry in range(start, end):
        tree.ReadEntry(entry)

        selectedJets = list(jetSelection(branchJet, ptcut=25, etacut=2.5))
        if len(selectedJets) < 6:
            continue
        bjets = list(bjetfilter(selectedJets))
        ljets = list(ljetfilter(selectedJets)) 

        electrons = list(branchElectron)
        muons = list(branchMuon)
        leptons = electrons + muons
        MET = branchMET[0]

        len_list = [len(leptons), len(selectedJets), len(ljets), len(bjets)]
        category = decay_channels(len_list)

        # bjet length between 2 and 4 produce most accurate reconstructions.
        if 2 <= len_list[3] <= 4:
            execute(category, HW, Ht, Hphi, selectedJets, leptons, MET)
        
      
def main():
    
    HW = ROOT.TH1F("W Boson", "W Reconstruction; Mass(GeV); Abundance", 80, 0, 200)
    Ht = ROOT.TH1F("Top Quark", "Top Reconstruction; Mass(GeV); Abundance", 80, 0, 400)
    Hphi = ROOT.TH1F("Phide Boson", "Phide Reconstruction; Mass(GeV); Abundance", 120, 0, 600) #120, 0, 600) 100, 0, 1000)
    
    chain = ROOT.TChain("Delphes")
    chain.Add('tag_1_delphes_events100.root')
    tree = ROOT.ExRootTreeReader(chain)
    
    branchJet = tree.UseBranch("Jet")
    branchElectron = tree.UseBranch("Electron")
    branchMuon = tree.UseBranch("Muon")
    branchMET = tree.UseBranch("MissingET")

    analyze(0, tree.GetEntries(), tree, branchJet, branchElectron, branchMuon, branchMET, HW, Ht, Hphi)
            
    with ROOT.TFile("phide100_nobias.root", "recreate") as file:
        file.WriteObject(HW, "W")
        file.WriteObject(Ht, "top")
        file.WriteObject(Hphi, "phide")
    
    c = ROOT.TCanvas()
    c.Divide(1,2)
    c.cd(1), Ht.Draw()
    c.cd(2), HW.Draw()
    c.cd(3), Hphi.Draw()
    c.SaveAs("phide100_nobias.pdf")

    end_time = time()
    print(f"Runtime: {end_time - start_time} seconds")

if __name__ == '__main__':
    main()

