#!/usr/bin/env python
# -*- coding: utf-8 -*-

import tdrstyle

runs = ['A', 'B', 'C', 'D']#, 'ABC', 'ABCD']

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
root_files_date = '2019-09-03'
root_files_runs = {}
for run in runs:
    PhotonJetPlots_dir_str = 'PhotonJetPlots_DATA_only_L2Res_0_3_.*_{run}_only_L2Res_vs_MC_only_L2Res_0_3_.*_{run}_only_L2Res_PFlowAK4chs_LUMI_vs_pt'.format(run=run)
    root_files_runs[run] = '{}/{}/plots.root'.format(root_files_base_dir, get_most_recent(root_files_base_dir, cut_str = PhotonJetPlots_dir_str))
for key, file in root_files_runs.iteritems():
    root_files_runs[key] = TFile(file, 'READ')

lumis_runs = {
    'A' : 13654.355526985,
    'B' : 7057.825158567,
    'C' : 6894.770971269,
    'D' : 31066.589629726,
}
for run in ['ABC', 'ABCD']:
    lumis_runs[run] = 0
    for r in run:
        lumis_runs[run] += lumis_runs[r]

lumi_total = lumis_runs['ABCD']

xsec_pb = {
    'MC': 283000.0,
}

etas = ['00_08', '08_13', '13_19', '19_25', '25_30', '30_32', '32_52', '00_13', '00_03', '03_05', '05_08', '08_10', '10_13', '13_15', '15_17', '17_19', '19_22', '22_23', '23_25', '25_27', '27_29', '29_30', '30_31', '31_35', '35_38', '38_52']
etas = ['00_13']

GraphErrors = {}
GraphHistos = {}
for key in root_files_runs:
    GraphErrors[key] = {}
    GraphHistos[key] = {}
    for eta in etas:
        GraphErrors[key][eta] = root_files_runs[key].Get('resp_MPFchs_DATA_a30_eta{}'.format(eta))
GraphErrors['MC'] = {}
GraphHistos['MC'] = {}
for eta in etas:
    GraphErrors['MC'][eta] = root_files_runs['A'].Get('resp_MPFchs_MC_a30_eta{}'.format(eta))

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

os.system("mkdir -p {}/MPF_per_runs_plots/".format(root_files_base_dir))

colors = {
    'A' : 1,
    'B' : 2,
    'C' : 4,
    'D' : 6,
    'MC' : 3
}
FakeGraphErrors = GraphErrors['MC'][etas[0]].Clone()

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

for eta in etas:
    can, pad, padr = buildCanvas(eta)
    pad.cd()
    FakeGraphErrors = GraphErrors['MC'][eta].Clone()
    FakeGraphErrors.SetPointError(FakeGraphErrors.GetN()-1, Double(4000), Double(0))
    y_err = 0
    for run in ['MC']+runs:
        if GraphErrors[run][eta].GetErrorY(0) > y_err:
            y_err = GraphErrors[run][eta].GetErrorY(0)
    if y_err > 0.03:
        rm_first_point = True
    else:
        rm_first_point = False
    y_err = 0
    for run in ['MC']+runs:
        if GraphErrors[run][eta].GetErrorY(1) > y_err:
            y_err = GraphErrors[run][eta].GetErrorY(1)
    if y_err > 0.03:
        rm_second_point = True
    else:
        rm_second_point = False

    FakeGraphErrors.SetMarkerSize(0)
    FakeGraphErrors.SetLineWidth(0)
    FakeGraphErrors.SetTitle("")
    if rm_first_point:
        remove_first_point(FakeGraphErrors)
    if rm_second_point:
        remove_second_point(FakeGraphErrors)
    FakeGraphErrors.Draw("")

    legend = TLegend(.75,.6,.925,.9)
    legend.SetLineColor(0)
    legend.Draw()

    GraphErrors['MC'][eta].SetMarkerColor(colors['MC'])
    GraphErrors['MC'][eta].SetMarkerStyle(24)
    GraphErrors['MC'][eta].SetMarkerSize(1)
    GraphErrors['MC'][eta].SetLineWidth(1)
    if rm_first_point:
        remove_first_point(GraphErrors['MC'][eta])
    if rm_second_point:
        remove_second_point(GraphErrors['MC'][eta])
    GraphErrors['MC'][eta].Draw("P same")
    add_legend_entry(GraphErrors['MC'][eta], 'MC')
    for run in runs:
        GraphErrors[run][eta].SetMarkerColor(colors[run])
        GraphErrors[run][eta].SetMarkerStyle(20)
        GraphErrors[run][eta].SetMarkerSize(1)
        GraphErrors[run][eta].SetLineWidth(1)
        if rm_first_point:
            remove_first_point(GraphErrors[run][eta])
        if rm_second_point:
            remove_second_point(GraphErrors[run][eta])
        GraphErrors[run][eta].Draw("P same")
        add_legend_entry(GraphErrors[run][eta], 'Run {}'.format(run))

    legend.Draw()

    Xaxis = FakeGraphErrors.GetXaxis()
    Xaxis.SetRangeUser(40,2000)
    pad.SetLogx()
    Xaxis.SetMoreLogLabels()
    Xaxis.SetNoExponent()
    Xaxis.SetLabelSize(0)

    Yaxis = FakeGraphErrors.GetYaxis()
    Yaxis.SetRangeUser(.7,1.4)
    Yaxis.SetTitle('MPF')
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
    pave3.AddText("2018 (13 TeV)")
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

    same = ""
    TgraphratioMC = GraphErrors['MC'][eta]
    Tgraphratios = []
    for run in runs:
        Tgraphratios.append(GraphErrors[run][eta].Clone())
        for k in range(Tgraphratios[-1].GetN()):
            x_value = Double(1)
            y_value = Double(1)
            x_valueMC = Double(1)
            y_valueMC = Double(1)
            Tgraphratios[-1].GetPoint(k, x_value, y_value)
            TgraphratioMC.GetPoint(k, x_valueMC, y_valueMC)
            if y_valueMC > 0:
                Tgraphratios[-1].SetPoint(k, x_value, y_value/y_valueMC)
                Tgraphratios[-1].SetPointError(k, Double(0),
                                      Double(Tgraphratios[-1].GetErrorY(k)/y_valueMC))
            else:
                Tgraphratios[-1].SetPoint(k, x_value, Double(1))
                Tgraphratios[-1].SetPointError(k, Double(0),
                                               Double(0))
                
        if same is not "P same":
            ratioFakeGraphErrors = Tgraphratios[-1].Clone()
            for k in range(ratioFakeGraphErrors.GetN()):
                x_value = Double(1)
                y_value = Double(1)
                ratioFakeGraphErrors.GetPoint(k, x_value, y_value)
                ratioFakeGraphErrors.SetPoint(k, x_value, Double(1))
                ratioFakeGraphErrors.SetPointError(k, Double(0), Double(0))
            ratioFakeGraphErrors.SetPointError(k, Double(4000), Double(0))
            ratioFakeGraphErrors.SetMarkerSize(0)
            ratioFakeGraphErrors.SetLineWidth(1)
            ratioFakeGraphErrors.SetTitle("")
            ratioFakeGraphErrors.Draw("")
            same = "P same"
        Tgraphratios[-1].Draw(same)

    ratioXaxis = ratioFakeGraphErrors.GetXaxis()
    ratioXaxis.SetRangeUser(40,2000)
    ratioXaxis.SetTitle('$p_T^\gamma$ \\text{(GeV)}')
    padr.SetLogx()
    ratioXaxis.SetMoreLogLabels()
    ratioXaxis.SetNoExponent()
    ratioXaxis.SetLabelSize(.075)
    ratioXaxis.SetTitleSize(.075)

    ratioYaxis = ratioFakeGraphErrors.GetYaxis()
    ratioYaxis.SetRangeUser(.89,1.11)
    ratioYaxis.SetTitle('Data/MC')
    ratioYaxis.SetTitleSize(.075)
    ratioYaxis.SetTitleOffset(.55)
    ratioYaxis.SetLabelSize(.075)
    ratioYaxis.SetNdivisions(5)

    padr.Update()
    for ext in ['pdf', 'png', 'C', 'root', 'tex']:
        can.SaveAs("{}/MPF_per_runs_plots/resp_MPF_eta{}_all_runs.{}".format(root_files_base_dir,eta,ext))
