from ROOT import *

#f = TFile("PhotonJet_2ndLevel_SinglePhoton_RunB2015.root")
f = TFile("../tuples/GJET_Pythia/PhotonJet_2ndLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8__25ns_ReReco_2016-05-31.root")

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
