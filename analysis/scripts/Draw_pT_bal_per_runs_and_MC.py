#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
user = os.popen("echo $USER").readlines()[0]
CMSSW_BASE = os.popen("echo $CMSSW_BASE").readlines()[0][:-1]

def get_most_recent(directory, cut_str = None):
    if cut_str is None:
        return os.popen('ls -t {}'.format(directory)).read().split('\n')[0]
    else:
        return os.popen('ls -t {} | grep {}'.format(directory, cut_str)).read().split('\n')[0]

from ROOT import TH1F, TFile, TCanvas, Double

user = "ltortero"
root_files_base_dir = "{}/src/JetMETCorrections/GammaJetFilter".format(CMSSW_BASE)
root_files_date = '2019-09-03'
root_files_runs = {}
for run in ['A', 'B', 'C', 'D', 'ABC', 'ABCD']:
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

GraphErrors = {}
for key in root_files_runs:
    GraphErrors[key] = {}
    for eta in etas:
        GraphErrors[key][eta] = root_files_runs[key].Get('resp_PtBalchs_DATA_a30_eta{}'.format(eta))
GraphErrors['MC'] = {}
for eta in etas:
    GraphErrors['MC'][eta] = root_files_runs['A'].Get('resp_PtBalchs_MC_a30_eta{}'.format(eta))

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
    pad.SetBottomMargin(0.03)
    pad.SetRightMargin(0.05)
    
    padr.SetBottomMargin(0.35)
    padr.SetLeftMargin(0.16)
    padr.SetRightMargin(0.05)
        
    can.cd()
    import locale; locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
    can.Draw()
    pad.Draw()
    padr.Draw()

    return can, pad, padr

os.system("mkdir -p {}/pT_bal__per_runs_plots/".format(root_files_base_dir))

colors = {
    'A' : 1,
    'B' : 2,
    'C' : 4,
    'D' : 6,
    'MC' : 3
}

for eta in etas:
    can, pad, padr = buildCanvas(eta)
    pad.cd()
    GraphErrors['MC'][eta].SetMarkerColor(colors['MC'])
    GraphErrors['MC'][eta].SetMarkerStyle(24)
    GraphErrors['MC'][eta].SetLineWidth(1)
    GraphErrors['MC'][eta].Draw("")
    for run in ['A', 'B', 'C', 'D']:
        GraphErrors[run][eta].SetMarkerColor(colors[run])
        GraphErrors[run][eta].SetMarkerStyle(20)
        GraphErrors[run][eta].SetLineWidth(1)
        GraphErrors[run][eta].Draw("SAME")

    Xaxis = GraphErrors['MC'][eta].GetXaxis()
    Xaxis.SetRangeUser(40,2000)
    pad.SetLogx()

    Yaxis = GraphErrors['MC'][eta].GetYaxis()
    Yaxis.SetRangeUser(.8,1.2)
    Yaxis.SetTitle('p_T balance')

    pad.Update()
    padr.cd()

    is_same = ""
    TgraphratioMC = GraphErrors['MC'][eta]
    for run in ['A', 'B', 'C', 'D']:
        Tgraphratio = GraphErrors[run][eta].Clone()
        for k in range(Tgraphratio.GetN()):
            x_value = Double(1)
            y_value = Double(1)
            x_valueMC = Double(1)
            y_valueMC = Double(1)
            Tgraphratio.GetPoint(k, x_value, y_value)
            TgraphratioMC.GetPoint(k, x_valueMC, y_valueMC)
            Tgraphratio.SetPoint(k, x_value, y_value/y_valueMC)
            Tgraphratio.SetPointError(k, Double(0),
                                      Double(Tgraphratio.GetErrorY(k)/y_valueMC))
        Tgraphratio.Draw(is_same)
        is_same = "SAME"

    ratioXaxis = Tgraphratio.GetXaxis()
    ratioXaxis.SetRangeUser(40,2000)
    ratioXaxis.SetTitle('p_T^#gamma (GeV)')
    padr.SetLogx()
    ratioXaxis.SetMoreLogLabels()
    ratioXaxis.SetNoExponent()

    ratioYaxis = Tgraphratio.GetYaxis()
    ratioYaxis.SetRangeUser(.95,1.05)
    ratioYaxis.SetTitle('Data/MC')

    padr.Update()

    import pdb; pdb.set_trace()
# for process in shapes:
#     print process
#     integrals_infos[process] = {}
#     Nb_events_infos[process] = {}
#     integrals_max_ratios[process] = {}
#     integrals_min_ratios[process] = {}
#     integrals_max_ratios_err[process] = {}
#     integrals_min_ratios_err[process] = {}
#     for syst in shapes[process]:
#         integrals_infos[process][syst] = {}
#         Nb_events_infos[process][syst] = {}
#         integrals_max_ratios[process][syst] = 1
#         integrals_min_ratios[process][syst] = 1
#         integrals_max_ratios_err[process][syst] = 0
#         integrals_min_ratios_err[process][syst] = 0
#         Nb_bins_for_infos += 1
#         histo_nomi = histdir.Get(process)
#         try:
#             histo_up = histdir.Get(shapes[process][syst]['Up'])
#         except(KeyError):
#             print "No Up for {} {}".format(process, syst)
#             histo_up = copy.copy(histo_nomi)
#         try:
#             histo_down = histdir.Get(shapes[process][syst]['Down'])
#         except(KeyError):
#             print "No Down for {} {}".format(process, syst)
#             histo_down = copy.copy(histo_nomi)

#         histo_nomi.SetFillColor(0)
#         histo_up.SetFillColor(0)
#         histo_down.SetFillColor(0)
#         histo_nomi.SetMarkerStyle(2)
#         histo_up.SetMarkerStyle(26)
#         histo_down.SetMarkerStyle(24)
#         histo_nomi.SetMarkerColor(1)
#         histo_nomi.SetLineColor(1)
#         histo_up.SetMarkerColor(styles.EWKcol)
#         histo_up.SetLineColor(styles.EWKcol)
#         histo_down.SetMarkerColor(styles.dycol)
#         histo_down.SetLineColor(styles.dycol)
        
#         can, pad, padr = buildCanvas(process, syst=syst)
#         pad.cd()

#         histo_nomi.SetTitle('{} {}'.format(process, syst))
#         histo_nomi.SetStats(0)
#         histo_nomi.Draw('')
#         histo_up.SetTitle('')
#         histo_down.SetTitle('')
#         histo_up.Draw("same")
#         histo_down.Draw("same")
#         histo_nomi.Draw("same")

#         Xaxis = histo_nomi.GetXaxis()
#         Yaxis = histo_nomi.GetYaxis()
#         if 'btag' in options.category and not 'nobtag' in options.category:
#             Xaxis.SetRangeUser(20,4000)
#         else:
#             Xaxis.SetRangeUser(10,4000)
#         Yaxis.SetTitle('N events (weighted)')
#         pad.Update()
        
#         integrals_infos[process][syst]['Up'] = round(histo_up.Integral(),2)
#         Nb_events_infos[process][syst]['Up'] = histo_up.GetEntries()
        
#         integrals_infos[process][syst]['Down'] = round(histo_down.Integral(),2)
#         Nb_events_infos[process][syst]['Down'] = histo_down.GetEntries()
        
#         integrals_infos[process][syst]['Nomi'] = round(histo_nomi.Integral(),2)
#         Nb_events_infos[process][syst]['Nomi'] = histo_nomi.GetEntries()
        
#         pave = ROOT.TPaveText(.75, .4, .925, .65, 'ndc')
#         for typ in ['Up', 'Nomi', 'Down']:
#             pave.AddText("int {}: {}".format(typ, integrals_infos[process][syst][typ]))
#         for typ in ['Up', 'Nomi', 'Down']:
#             pave.AddText("ent {}: {}".format(typ, Nb_events_infos[process][syst][typ]))
#         pave.SetTextSizePixels(15)
#         pave.SetTextAlign(11)
#         pave.SetBorderSize(0)
#         pave.SetFillColor(0)
#         pave.SetFillStyle(0)
#         pave.Draw()

#         legend = ROOT.TLegend(.75,.7,.925,.9)
#         legend.AddEntry(histo_up.GetName(), 'Up')
#         legend.AddEntry(histo_nomi.GetName(), 'Nominal')
#         legend.AddEntry(histo_down.GetName(), 'Down')
#         legend.Draw()
        
#         padr.cd()

#         rhisto_down = copy.copy(histo_down)
#         rhisto_up = copy.copy(histo_up)
#         rhisto_down.Divide(histo_nomi)
#         rhisto_up.Divide(histo_nomi)

#         for b in range(histo_nomi.GetNbinsX()):
#             rvalue_up = rhisto_up.GetBinContent(b+1)
#             rvalue_down = rhisto_down.GetBinContent(b+1)
#             if rvalue_up > integrals_max_ratios[process][syst]:
#                 integrals_max_ratios[process][syst] = rvalue_up
#                 integrals_max_ratios_err[process][syst] = rhisto_up.GetBinError(b+1)
#             if rvalue_down > integrals_max_ratios[process][syst]:
#                 integrals_max_ratios[process][syst] = rvalue_down
#                 integrals_max_ratios_err[process][syst] = rhisto_down.GetBinError(b+1)
#             if rvalue_up < integrals_min_ratios[process][syst] and rvalue_up > 0:
#                 integrals_min_ratios[process][syst] = rvalue_up
#                 integrals_min_ratios_err[process][syst] = rhisto_up.GetBinError(b+1)
#             if rvalue_down < integrals_min_ratios[process][syst] and rvalue_down > 0:
#                 integrals_min_ratios[process][syst] = rvalue_down
#                 integrals_min_ratios_err[process][syst] = rhisto_down.GetBinError(b+1)

#         rhisto_down.SetStats(0)
#         rhisto_down.Draw()
#         rhisto_up.Draw("same")

#         ratioXaxis = rhisto_down.GetXaxis()
#         ratioYaxis = rhisto_down.GetYaxis()
#         if 'btag' in options.category and not 'nobtag' in options.category:
#             ratioXaxis.SetRangeUser(20,4000)
#         else:
#             ratioXaxis.SetRangeUser(10,4000)
#         ratioYaxis.SetRangeUser(.5,1.5)
        
#         ratioXaxis.SetTitle('mT tot (GeV)')
#         ratioXaxis.SetTitleSize(.1)
        
#         ratioYaxis.SetTitle('Ratio / Nominal')
#         ratioYaxis.SetTitleSize(.1)
#         ratioYaxis.SetTitleOffset(0.5)
        
#         ratioXaxis.SetMoreLogLabels()
#         ratioXaxis.SetNoExponent()

#         pad.SetLogx()
#         padr.SetLogx()
#         padr.Update()

#         # legend = ROOT.TLegend(.75,.9,1,1)
#         # legend.AddEntry('up_{}'.format(infos), 'Up/Nominal')
#         # legend.AddEntry('do_{}'.format(infos), 'Down/Nominal')
#         # legend.Draw()
        
#         syst = syst.replace('CMS_htt_{}_{}_{}_13TeV__'.format(options.channel, options.channel, options.category),'')
#         can.SaveAs("./shapes_ctrl_plots/{}/{}/{}_{}.png".format(options.channel, options.category, process, syst))

# ratios = {'integrals':integrals_infos,'entries':Nb_events_infos}
# for infos in ratios:
#     can = ROOT.TCanvas('canratios_{}'.format(infos), '', 2000, 1500)
#     pad = can.GetPad(0)
    
#     # Set Pad sizes
#     pad.SetPad(0.0, 0.32, 1., 1.0)
    
#     pad.SetLeftMargin(.05)
#     pad.SetBottomMargin(.5)
#     pad.SetRightMargin(.01)
        
#     can.cd()
#     import locale; locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
#     can.Draw()
#     pad.Draw('same')
    
#     pad.cd()
    
#     binning = Nb_bins_for_infos
#     histo_up = ROOT.TH1F('up_{}'.format(infos), 'Ratios of {}'.format(infos), *(binning,0,binning))
#     histo_down = ROOT.TH1F('do_{}'.format(infos), 'Ratios of {}'.format(infos), *(binning,0,binning))
    
#     histo_up.SetMarkerStyle(26)
#     histo_down.SetMarkerStyle(24)
#     histo_up.SetMarkerColor(styles.EWKcol)
#     histo_up.SetLineColor(styles.EWKcol)
#     histo_down.SetMarkerColor(styles.dycol)
#     histo_down.SetLineColor(styles.dycol)
    
#     bin = 0
#     for process in ratios[infos]:
#         for syst in ratios[infos][process]:
#             bin += 1
#             proc_syst = ratios[infos][process][syst]
#             histo_up.SetBinContent(bin, proc_syst['Up']/proc_syst['Nomi'])
#             histo_up.SetBinError(bin, proc_syst['Up']**(.5)/proc_syst['Nomi'])
#             histo_up.GetXaxis().SetBinLabel(bin,'{} {}'.format(syst, process))
#             histo_down.SetBinContent(bin, proc_syst['Down']/proc_syst['Nomi'])
#             histo_down.SetBinError(bin, proc_syst['Down']**(.5)/proc_syst['Nomi'])
#     histo_up.GetXaxis().SetLabelSize(0.025)
#     histo_up.SetStats(0)
#     histo_up.Draw()
#     histo_down.Draw('same')

#     legend = ROOT.TLegend(.75,.9,1,1)
#     legend.AddEntry('up_{}'.format(infos), 'Up/Nominal')
#     legend.AddEntry('do_{}'.format(infos), 'Down/Nominal')
#     legend.Draw()
    
#     can.SaveAs("./shapes_ctrl_plots/{}/{}/ratios_{}.png".format(options.channel, options.category, infos))

    
# can = ROOT.TCanvas('canratios', '', 2000, 1500)
# pad = can.GetPad(0)
    
# # Set Pad sizes
# pad.SetPad(0.0, 0.32, 1., 1.0)

# pad.SetLeftMargin(.05)
# pad.SetBottomMargin(.5)
# pad.SetRightMargin(.01)
        
# can.cd()
# import locale; locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
# can.Draw()
# pad.Draw('same')

# pad.cd()

# binning = Nb_bins_for_infos
# histo_max = ROOT.TH1F('max_vars'.format(infos), 'Relative variations', *(binning,0,binning))
# histo_min = ROOT.TH1F('min_vars'.format(infos), '', *(binning,0,binning))

# bin = 0
# for process in integrals_max_ratios:
#     for syst in integrals_max_ratios[process]:
#         bin += 1
#         histo_max.SetBinContent(bin, integrals_max_ratios[process][syst])
#         histo_max.SetBinError(bin, integrals_max_ratios_err[process][syst])
#         histo_min.SetBinContent(bin, integrals_min_ratios[process][syst])
#         histo_min.SetBinError(bin, integrals_min_ratios_err[process][syst])
#         histo_max.GetXaxis().SetBinLabel(bin,'{} {}'.format(syst, process))
#     histo_max.GetXaxis().SetLabelSize(0.025)
#     histo_max.GetYaxis().SetRangeUser(.5,1.5)
#     histo_max.SetStats(0)
#     histo_max.Draw()
#     histo_min.Draw('same')

# # legend = ROOT.TLegend(.75,.9,1,1)
# # legend.AddEntry('up_{}'.format(infos), 'Up/Nominal')
# # legend.AddEntry('do_{}'.format(infos), 'Down/Nominal')
# # legend.Draw()

# can.SaveAs("./shapes_ctrl_plots/{}/{}/relative_variations.png".format(options.channel, options.category))

    
# data_file.Close()
