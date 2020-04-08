#!/usr/bin/env python
# -*- coding: utf-8 -*-

import tdrstyle

runs = ['B', 'C', 'D', 'E', 'F']
runs+= ['BC', 'DE', 'BCDEF']

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
root_files_base_dir = "{}/src/JetMETCorrections/GammaJetFilter".format(CMSSW_BASE)

L1Offset_approaches = ["ComplexL1", "SimpleL1"]
JERCs = ['wo_L2Res', 'only_L2Res', 'L2L3Res', 'JER']

root_files = {}
allowed = []

for L1Offset_approach in L1Offset_approaches:
    root_files[L1Offset_approach] = {}
    for JERC in JERCs:
        root_files[L1Offset_approach][JERC] = {}
        for run in runs:
            root_files[L1Offset_approach][JERC][run] = None

            PhotonJetPlots_dir_str = 'PhotonJetPlots_DATA_{JERC}_0_3_.*_{run}_{JERC}_vs_MC_{JERC}_0_3_.*_{run}_{JERC}_PFlowAK4chs_LUMI_vs_pt'.format(run=run, JERC = '_'.join([L1Offset_approach,JERC]))
            file = '{}/{}/plots.root'.format(root_files_base_dir, get_most_recent(root_files_base_dir, cut_str =
PhotonJetPlots_dir_str))
            if file != "{}//plots.root".format(root_files_base_dir) :
                #print(file)
                root_files[L1Offset_approach][JERC][run] = file
                allowed.append("/".join([L1Offset_approach,JERC,run]))

lumis_runs = {
    'B' : 4793.961426839,
    'C' : 9631.214820913,
    'D' : 4247.682053046,
    'E' : 9313.642401775,
    'F' : 13539.378417564,
}

lumi_total = 0
for run in lumis_runs.keys():
    lumi_total += lumis_runs[run]

for run in runs:
    if run not in lumis_runs.keys():
        lumis_runs[run] = 0
        for r in run:
            lumis_runs[run] += lumis_runs[r]

xsec_pb = {
    'HT-40To100' : 18700.0,
    'HT-100To200' : 8640.0,
    'HT-200To400' : 2185.0,
    'HT-400To600' : 259.9,
    'HT-600ToInf': 85.31,
}

etas = ['00_08', '08_13', '13_19', '19_25', '25_30', '30_32', '32_52', '00_13', '00_03', '03_05', '05_08', '08_10', '10_13', '13_15', '15_17', '17_19', '19_22', '22_23', '23_25', '25_27', '27_29', '29_30', '30_31', '31_35', '35_38', '38_52']

responses = {
    "MPF" : "MPF",
    "BAL" : "PtBal",
}

def buildCanvas(name):
    can = TCanvas('can'+name, '', 800, 800)
    can.Divide(1, 2, 0.0, 0.0)
    
    pad = can.GetPad(1)
    padr = can.GetPad(2)
    
    # Set Pad sizes
    pad.SetPad(0.0, 0.32, 1., 1.0)
    padr.SetPad(0.0, 0.00, 1., 0.34)
    
    pad.SetTopMargin(0.08)
    pad.SetLeftMargin(0.16)
    pad.SetBottomMargin(0.05)
    pad.SetRightMargin(0.05)
    
    padr.SetBottomMargin(0.25)
    padr.SetLeftMargin(0.16)
    padr.SetRightMargin(0.05)
        
    can.cd()
    import locale; locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
    can.Draw()
    pad.Draw()
    padr.Draw()

    return can, pad, padr

outdir = "{}/Mikkos_comparison_plots/".format(root_files_base_dir)
os.system("mkdir -p {}".format(outdir))

colors = {
    'SimpleL1' : 1,
    'ComplexL1' : 2,
}

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


# GraphErrors = {}
# GraphHistos = {}
# for key in root_files_runs:
#     GraphErrors[key] = {}
#     GraphHistos[key] = {}
#     for eta in etas:
#         GraphErrors[key][eta] = root_files_runs[key].Get('resp_PtBalchs_DATA_a30_eta{}'.format(eta))
# GraphErrors['MC'] = {}
# GraphHistos['MC'] = {}
# for eta in etas:
#     GraphErrors['MC'][eta] = root_files_runs['A'].Get('resp_PtBalchs_MC_a30_eta{}'.format(eta))

for run in runs:
    for JERC in JERCs:
        proceed = all(["/".join([L1Offset_approach,JERC,run]) in allowed for L1Offset_approach in L1Offset_approaches])
        if not proceed:
            continue
        root_file = {}
        for L1Offset_approach in L1Offset_approaches:
            root_file[L1Offset_approach] = TFile(root_files[L1Offset_approach][JERC][run], 'READ')
        for eta in etas:
            if eta != '00_13':
                continue
            for response in responses:
                GraphErrorData = {}
                GraphErrorMC = {}
                for L1Offset_approach in L1Offset_approaches:
                    GraphErrorData[L1Offset_approach] = root_file[L1Offset_approach].Get('resp_{}chs_DATA_a30_eta{}'.format(responses[response],eta))
                    GraphErrorMC[L1Offset_approach]  =  root_file[L1Offset_approach].Get('resp_{}chs_MC_a30_eta{}'.format(responses[response],eta))

                for plottype in ["data", "mc", "data_over_mc"]:
                    can, pad, padr = buildCanvas("_".join([JERC,run,eta,response,plottype]))
                    pad.cd()

                    if plottype == "data":
                        graphs = GraphErrorData
                        y_title_str = "2017UL {}".format(run)
                    if plottype == "mc":
                        graphs = GraphErrorMC
                        y_title_str = "MC"
                    if plottype == "data_over_mc":
                        graphs = {}
                        y_title_str = "2017UL {} / MC".format(run)
                        for key in GraphErrorData.keys():
                            graphs[key] = GraphErrorData[key].Clone()
                            for k in range(graphs[key].GetN()):
                                x_value = Double(1)
                                y_value = Double(1)
                                x_valueMC = Double(1)
                                y_valueMC = Double(1)
                                GraphErrorData[key].GetPoint(k, x_value, y_value)
                                GraphErrorMC[key].GetPoint(k, x_valueMC, y_valueMC)
                                if y_valueMC > 0:
                                    graphs[key].SetPoint(k, x_value, y_value/y_valueMC)
                                    graphs[key].SetPointError(k, Double(0),Double(GraphErrorData[key].GetErrorY(k)/y_valueMC))
                                else:
                                    graphs[key].SetPoint(k, x_value, Double(1))
                                    graphs[key].SetPointError(k, Double(0),Double(0))
                    print(run, JERC, eta, response, plottype)
                    graphs["SimpleL1"].SetMarkerColor(colors['SimpleL1'])
                    graphs["SimpleL1"].SetMarkerStyle(24)
                    graphs["SimpleL1"].SetMarkerSize(1)
                    graphs["SimpleL1"].SetLineWidth(1)
                    graphs["SimpleL1"].Draw("")
                    graphs["ComplexL1"].SetMarkerColor(colors['ComplexL1'])
                    graphs["ComplexL1"].SetMarkerStyle(24)
                    graphs["ComplexL1"].SetMarkerSize(1)
                    graphs["ComplexL1"].SetLineWidth(1)
                    graphs["ComplexL1"].Draw("P same")

                    legend = TLegend(.75,.6,.925,.9)
                    legend.SetLineColor(0)
                    legend.Draw()

                    for L1Offset_approach in L1Offset_approaches:
                        add_legend_entry(graphs[L1Offset_approach], L1Offset_approach)

                    legend.Draw()
                        
                    Xaxis = graphs["SimpleL1"].GetXaxis()
                    Xaxis.SetRangeUser(10,4000)
                    pad.SetLogx()
                    Xaxis.SetMoreLogLabels()
                    Xaxis.SetNoExponent()
                    Xaxis.SetLabelSize(0)

                    Yaxis = graphs["SimpleL1"].GetYaxis()
                    Yaxis.SetRangeUser(.7,1.4)
                    Yaxis.SetTitle("Resp. {}, {}".format(response, y_title_str))
                    Yaxis.SetTitleSize(.05)

                    pave = TPaveText(.4, .6, .925-.35, .75, 'ndc')
                    pave.SetFillColor(0)
                    pave.SetFillStyle(0)
                    pave.SetBorderSize(0)
                    pave.AddText("$\\alpha < 0.3$")
                    pave.AddText("${} < |\eta| < {}$".format(int(eta.split('_')[0])*1./10, int(eta.split('_')[1])*1./10))
                    pave.Draw()

                    pave2 = TPaveText(0.15, 0.925, 0.5, 0.96, "ndc")
                    pave2.SetFillColor(0)
                    pave2.SetFillStyle(0)
                    pave2.SetBorderSize(0)
                    pave2.AddText("\gamma + \\text{jets}")
                    pave2.SetTextAlign(11)
                    pave2.SetTextSize(0.045)
                    pave2.Draw()

                    pave3 = TPaveText(0.5, 0.925, 0.96, 0.96, "ndc")
                    pave3.SetFillColor(0)
                    pave3.SetFillStyle(0)
                    pave3.SetBorderSize(0)
                    pave3.AddText("2017 UL (13 TeV)")
                    pave3.SetTextAlign(31)
                    pave3.SetTextSize(0.045)
                    pave3.Draw()

                    pave4 = TPaveText(0.175, 0.85, 0.5, 0.88, "brNDC")
                    pave4.SetFillColor(0)
                    pave4.SetFillStyle(0)
                    pave4.SetBorderSize(0)
                    pave4.AddText("CMS #it{Preliminary}")
                    pave4.SetTextAlign(11)
                    pave4.SetTextSize(0.045)
                    pave4.Draw()

                    pad.Update()
                    
                    padr.cd()
                    
                    rgraph = graphs["SimpleL1"].Clone()
                    for k in range(rgraph.GetN()):
                        x_value_a = Double(1)
                        y_value_a = Double(1)
                        x_value_b = Double(1)
                        y_value_b = Double(1)
                        graphs["SimpleL1"].GetPoint(k, x_value_a, y_value_a)
                        graphs["ComplexL1"].GetPoint(k, x_value_b, y_value_b)
                        if y_value_b > 0:
                            rgraph.SetPoint(k, x_value_a, y_value_a/y_value_b)
                            rgraph.SetPointError(k, Double(0),
                                                 Double(graphs["SimpleL1"].GetErrorY(k)/y_value_b))
                        else:
                            rgraph.SetPoint(k, x_value, Double(1))
                            rgraph.SetPointError(k, Double(0),
                                                 Double(0))
                
                    rgraph.SetMarkerSize(0)
                    rgraph.SetLineWidth(1)
                    rgraph.SetTitle("")
                    rgraph.Draw("")

                    ratioXaxis = rgraph.GetXaxis()
                    ratioXaxis.SetRangeUser(10,4000)
                    ratioXaxis.SetTitle('$p_T^\gamma$ \\text{(GeV)}')
                    padr.SetLogx()
                    ratioXaxis.SetMoreLogLabels()
                    ratioXaxis.SetNoExponent()
                    ratioXaxis.SetLabelSize(.075)
                    ratioXaxis.SetTitleSize(.075)
                    
                    ratioYaxis = rgraph.GetYaxis()
                    if plottype == "mc":
                        ratioYaxis.SetRangeUser(.99,1.01)                        
                    else:
                        ratioYaxis.SetRangeUser(.89,1.11)
                    ratioYaxis.SetTitle('Simple/Complex')
                    ratioYaxis.SetTitleSize(.075)
                    ratioYaxis.SetTitleOffset(.55)
                    ratioYaxis.SetLabelSize(.075)
                    ratioYaxis.SetNdivisions(5)

                    padr.Update()
                    for ext in ['pdf', 'png', 'C', 'root', 'tex']:
                        can.SaveAs("{}/resp_{}_{}_eta{}_run{}.{}".format(outdir, response, plottype, eta, run, ext))
        for L1Offset_approach in L1Offset_approaches:
            root_file[L1Offset_approach].Close()
