#!/usr/bin/env python
# -*- coding: utf-8 -*-

import tdrstyle

runs18 = ['A', 'B', 'C', 'D']
runs18+= ["ABC", "ABCD"]
runs17 = ['B', 'C', 'D', 'E', 'F']
runs17+= ['BC', 'DE', 'BCDEF']
HLTs = [33, 50, 75, 90, 120, 165, 200]

colors = {
    33 : 1,
    50 : 2,
    75 : 3,
    90 : 4,
    120 : 42,
    165 : 6,
    200 : 7,
}

markerstyles = {
    33 : 20,
    50 : 21,
    75 : 22,
    90 : 23,
    120 : 29,
    165 : 33,
    200 : 34,
}


import os
user = os.popen("echo $USER").readlines()[0]
CMSSW_BASE = os.popen("echo $CMSSW_BASE").readlines()[0][:-1]

def get_most_recent(directory, cut_str = None):
    if cut_str is None:
        return os.popen('ls -t {}'.format(directory)).read().split('\n')[0]
    else:
        return os.popen('ls -t {} | grep {}'.format(directory, cut_str)).read().split('\n')[0]

from ROOT import TFile, TCanvas, Double, TLegend, TH1F, TPaveText

user = "ltortero"
root_files_base_dir = "{}/src/JetMETCorrections/GammaJetFilter/analysis/PUReweighting/".format(CMSSW_BASE)

L1Offset_approaches = ["ComplexL1", "SimpleL1", "18"]
JERCs = ['wo_L2Res', 'only_L2Res', 'L2L3Res', 'JER']

root_files = {}
for L1Offset_approach in L1Offset_approaches:
    root_files[L1Offset_approach] = {}
    for JERC in JERCs:
        root_files[L1Offset_approach][JERC] = {}
        if L1Offset_approach == "18":
            for run in runs18:
                root_files[L1Offset_approach][JERC][run] = {}
                for HLT in HLTs:
                    root_files[L1Offset_approach][JERC][run][HLT] = "{root_files_base_dir}/pu_truth_data2018_100bins_HLTphoton{HLT}{run}_{JERC}.root".format(
                        root_files_base_dir=root_files_base_dir,
                        HLT=HLT,
                        run=run,
                        JERC=JERC)
        else:
            for run in runs17:
                root_files[L1Offset_approach][JERC][run] = {}
                for HLT in HLTs:
                    root_files[L1Offset_approach][JERC][run][HLT] = "{root_files_base_dir}/pu_truth_data2017UL_100bins_HLTphoton{HLT}{run}_{L1Offset_approach}_{JERC}.root".format(
                        root_files_base_dir=root_files_base_dir,
                        HLT=HLT,
                        run=run,
                        L1Offset_approach=L1Offset_approach,
                        JERC=JERC)

def buildCanvas(name):
    can = TCanvas('can'+name, '', 800, 800)
    can.Divide(1, 1, 0.0, 0.0)
    
    pad = can.GetPad(1)
    padr = None # can.GetPad(2)
    
    # Set Pad sizes
    #pad.SetPad(0.0, 0.32, 1., 1.0)
    #padr.SetPad(0.0, 0.00, 1., 0.34)
    
    pad.SetTopMargin(0.0)
    pad.SetLeftMargin(0.125)
    pad.SetBottomMargin(0.08)
    pad.SetRightMargin(0.0125)
    
    # padr.SetBottomMargin(0.25)
    # padr.SetLeftMargin(0.16)
    # padr.SetRightMargin(0.05)
        
    can.cd()
    import locale; locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
    can.Draw()
    pad.Draw()
    #padr.Draw()

    return can, pad, padr

outdir = "{}/PU_profiles_plots/".format(root_files_base_dir)
os.system("mkdir -p {}".format(outdir))

hfls = []
def add_legend_entry(GraphError, label):
    hfls.append(TH1F(label, label, 2, Double(1), Double(2)))
    hfls[-1].SetMarkerColor(GraphError.GetMarkerColor())
    hfls[-1].SetMarkerStyle(GraphError.GetMarkerStyle())
    hfls[-1].SetLineWidth(0)
    hfls[-1].Draw("same")
    legend.AddEntry(label, label)
    
def remove_point(tgrapherror, n=0):
    tgrapherror.SetPoint(n, Double(0), Double(0))
    tgrapherror.SetPointError(n, Double(0), Double(0))
    
def remove_first_point(tgrapherror):
    remove_point(tgrapherror, n=0)
    
def remove_second_point(tgrapherror):
    remove_point(tgrapherror, n=1)

for L1Offset_approach in L1Offset_approaches:
    if L1Offset_approach == "18":
        runs = runs18
        year = 2018
    else:
        runs = runs17
        year = 2017
    for run in runs:
        for JERC in JERCs:
            root_Tfiles = {}
            pileup_profiles = {}
            for HLT in HLTs:
                root_Tfiles[HLT] = TFile(root_files[L1Offset_approach][JERC][run][HLT], 'READ')
                pileup_profiles[HLT] = root_Tfiles[HLT].Get("pileup")
                ok = True
                try :
                    pileup_profiles[HLT].Scale(1/pileup_profiles[HLT].Integral())
                    pileup_profiles[HLT].SetMarkerColor(colors[HLT])
                    pileup_profiles[HLT].SetMarkerStyle(markerstyles[HLT])
                except AttributeError:
                    ok = False
            if not ok:
                continue
            can, pad, padr = buildCanvas("_".join([JERC,run]))
            pad.cd()
            same = ""
            for HLT in HLTs:
                pileup_profiles[HLT].Draw(same)
                same = "P same"

            y_title_str = "Probabilite"
            x_title_str = "Nombre d'interations d'empilement"
  
            legend = TLegend(.33,.125,.5,.6)
            legend.SetLineColor(0)
            legend.Draw()

            for HLT in HLTs:
                add_legend_entry(pileup_profiles[HLT], "HLT {}".format(str(int(HLT))))

            legend.Draw()
                        
            Xaxis = pileup_profiles[HLTs[0]].GetXaxis()
            Xaxis.SetTitle(x_title_str)
            #Xaxis.SetLabelSize(0)
            Xaxis.SetRangeUser(5,70)
        
            Yaxis = pileup_profiles[HLTs[0]].GetYaxis()
            Yaxis.SetTitleOffset(1.5)
            Yaxis.SetTitle(y_title_str)
            Yaxis.SetRangeUser(10**-5,1)
            pad.SetLogy()
            for ext in ['pdf', 'png', 'C', 'root', 'tex']:
                can.SaveAs("{outdir}/PU_HLT_profiles_run{year}{run}_{JERC}.{ext}".format(outdir=outdir, year=str(year), run=run, JERC=JERC, ext=ext))
            for tfile in root_Tfiles.values():
                tfile.Close()
       
