from ROOT import *

#f = TFile("PhotonJet_2ndLevel_SinglePhoton_RunB2015.root")
f = TFile("../tuples/QCD_MC/PhotonJet_2ndLevel_QCD_Pt-15toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8_25ns_ReReco_2015-12-03.root")

#f = TFile("PhotonJet_2ndLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_.root")
tree = f.Get("gammaJet/analysis")


for event in range(0,1):
    tree.GetEntry(event)
    size = tree.trigger_names.size()
    for i in range (0, size ):    
#        index = tree.triggerIndex(tree.trigger_names[trigger])
#        if tree.trigger_results[i] ==1 :
 #           passed = True
#        print event 
#        print i
        print tree.trigger_names[i]
#        print tree.trigger_results[i]
        print tree.trigger_prescale[i]
   #         print passed
