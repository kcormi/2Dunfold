# uncompyle6 version 3.9.0
# Python bytecode version base 2.7 (62211)
# Decompiled from: Python 3.7.1 (default, Dec 14 2018, 13:28:58) 
# [Clang 4.0.1 (tags/RELEASE_401/final)]
# Embedded file name: /work/jinw/CMSSW_10_2_15_patch2/src/cleanup_2Dunfold/plot_utils.py
# Compiled at: 2023-03-06 19:29:44
import numpy as np
from argparse import ArgumentParser
import os
from math import *
import ROOT as rt
import tdrstyle
from Plotting_cfg import *

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from unfold_utils import *

rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptStat(0)
tdrstyle.setTDRStyle()
rt.TH1.SetDefaultSumw2()
rt.TH2.SetDefaultSumw2()

colors = [rt.kBlue - 7, rt.kAzure + 5, rt.kCyan]
fillcolors = [38, 30, 42, 49, 14]
linecolors = [1, 2, 4, 7, 800, 880]

def DrawPlotStyle(HistList, PlotStyle, Color, Pad):
    Pad.cd()
    if PlotStyle == 'marker':
        for Hist in HistList:
            Hist.SetFillStyle(0)
            Hist.SetLineColor(Color)
            Hist.SetLineWidth(1)
            Hist.SetMarkerStyle(20)
            Hist.SetMarkerColor(Color)
            Hist.SetMarkerSize(0.5)
            Hist.DrawCopy('pesame')
            Pad.Update()

        LegendStyle = 'lep'
    if PlotStyle == 'cross':
        for Hist in HistList:
            Hist.SetFillStyle(0)
            Hist.SetLineColor(Color)
            Hist.SetLineWidth(1)
            Hist.SetMarkerStyle(20)
            Hist.SetMarkerColor(Color)
            Hist.SetMarkerSize(0)
            Hist.DrawCopy('pesame')
            Pad.Update()

        LegendStyle = 'lep'
    if PlotStyle == 'triangle':
        for Hist in HistList:
            Hist.SetFillStyle(0)
            Hist.SetLineColor(Color)
            Hist.SetLineWidth(1)
            Hist.SetMarkerStyle(26)
            Hist.SetMarkerColor(Color)
            Hist.SetMarkerSize(0.5)
            Hist.DrawCopy('pesame')
            Pad.Update()

        LegendStyle = 'lep'
    if PlotStyle == 'triangledown':
        for Hist in HistList:
            Hist.SetFillStyle(0)
            Hist.SetLineColor(Color)
            Hist.SetLineWidth(1)
            Hist.SetMarkerStyle(32)
            Hist.SetMarkerColor(Color)
            Hist.SetMarkerSize(0.5)
            Hist.DrawCopy('pesame')
            Pad.Update()

        LegendStyle = 'lep'
    if PlotStyle == 'hashed':
        for Hist in HistList:
            Hist.SetLineWidth(1)
            Hist.SetFillStyle(3004)
            Hist.SetFillColor(Color)
            Hist.SetLineColor(Color)
            Hist.SetMarkerSize(0)
            Hist.DrawCopy('e2same')
            Pad.Update()

        LegendStyle = 'f'
    if PlotStyle == 'hashedright':
        for Hist in HistList:
            Hist.SetLineWidth(1)
            Hist.SetFillStyle(3005)
            Hist.SetFillColor(Color)
            Hist.SetLineColor(Color)
            Hist.SetMarkerSize(0)
            Hist.DrawCopy('e2same')
            Pad.Update()

        LegendStyle = 'f'
    if PlotStyle == 'hist':
        for Hist in HistList:
            Hist.SetLineColor(Color)
            Hist.SetLineWidth(1)
            Hist.SetMarkerSize(0)
            Hist.SetFillStyle(0)
            Hist.DrawCopy('histsame')
            Pad.Update()

        LegendStyle = None
    if PlotStyle == 'histhashed':
        for Hist in HistList:
            Hist.SetLineWidth(1)
            Hist.SetFillStyle(3004)
            Hist.SetFillColor(Color)
            Hist.SetLineColor(Color)
            Hist.SetMarkerSize(0)
            Hist.DrawCopy('e2same')
            Pad.Update()
            Hist_step = Hist.Clone()
            Hist_step.SetLineWidth(1)
            Hist_step.SetLineColor(Color)
            Hist_step.SetFillStyle(0)
            Hist_step.SetMarkerSize(0)
            Hist_step.DrawCopy('histsame')
            Pad.Update()

        LegendStyle = 'f'
    if PlotStyle == 'filled':
        HistList_step = []
        for Hist in HistList:
            HistList_step.append(Hist.Clone())
            Hist.SetMarkerSize(0)
            Hist.SetFillStyle(1001)
            Hist.SetLineWidth(0)
            Hist.SetFillColorAlpha(Color, 0.75)
            Hist.DrawCopy('histsame')
            Pad.Update()

        for Hist_step in HistList_step:
            Hist_step.SetMarkerSize(0)
            Hist_step.SetFillStyle(1001)
            Hist_step.SetLineWidth(0)
            Hist_step.SetFillColorAlpha(Color, 0.95)
            Hist_step.DrawCopy('e2same')
            Pad.Update()

        LegendStyle = 'f'
    if PlotStyle == 'fillederror':
        HistList_step = []
        for Hist in HistList:
            HistList_step.append(Hist.Clone())

        for Hist_step in HistList_step:
            Hist_step.SetMarkerSize(0)
            Hist_step.SetFillStyle(0)
            Hist_step.SetLineWidth(1)
            Hist_step.SetLineColor(Color)
            Hist_step.DrawCopy('histsame')
            Pad.Update()

        for Hist in HistList:
            Hist.SetMarkerSize(0)
            Hist.SetFillStyle(1001)
            Hist.SetLineWidth(0)
            Hist.SetFillColorAlpha(Color, 0.95)
            Hist.DrawCopy('e2same')
            Pad.Update()

        LegendStyle = 'f'
    if PlotStyle == 'line':
        for Hist in HistList:
            Hist.SetLineColor(Color)
            Hist.SetLineWidth(1)
            Hist.SetMarkerStyle(21)
            Hist.SetMarkerColor(Color)
            Hist.SetMarkerSize(0.5)
            Hist.SetFillStyle(0)
            Hist.DrawCopy('lpsame')
            Pad.Update()

        LegendStyle = 'lp'
    Pad.Update()
    return LegendStyle


def plot_flat_hists(hist_ref, list_hist_compare, legend_ref, list_legend_compare, title, is_logY, do_ratio, output_path, hist_ref_stat=0, text_list=[], style_ref='marker', color_ref=rt.kBlack, list_style_compare=[], list_color_compare=[], labelY='Normalized Events/Bin Width', label_ratio='Data/MC'):
    W_merge = W_long
    if len(hist_ref.root_hists) > 3:
        W_merge = W_long * (float(len(hist_ref.root_hists)) / 3.0)
    if text_list == []:
        text_list = [
         ''] * len(hist_ref.root_hists)
    c = rt.TCanvas(f'c_{hist_ref.name}_{title}_{output_path}', f'c_{hist_ref.name}_{title}_{output_path}', int(W_merge), H_long)
    pad1 = rt.TPad('pad1', 'pad1', 0.0, 0.5 if do_ratio else 0.0, 1, 1)
    pad1.SetRightMargin(0.03)
    pad1.Draw()
    pad1.cd()
    pad1.SetLogy(is_logY)
    pad1.Update()
    xAxis = hist_ref.root_flat_hist.GetXaxis()
    xAxis.SetTitle(title)
    xAxis.SetTitleSize(FTS)
    xAxis.SetTitleOffset(1.05)
    if do_ratio:
        xAxis.SetLabelSize(0)
        xAxis.SetTitleSize(0)
        xAxis.SetTitleOffset(0)
    else:
        xAxis.SetLabelSize(FLS)
    yAxis = hist_ref.root_flat_hist.GetYaxis()
    yAxis.SetNdivisions(6, 5, 0)
    yAxis.SetLabelSize(FLS)
    yAxis.SetTitleSize(FTS)
    yAxis.SetMaxDigits(3)
    yAxis.SetTitle(labelY)
    yAxis.SetTitleOffset(0.5 + 0.8 / len(hist_ref.root_hists))
    yAxis.SetTickLength(yAxis.GetTickLength() / float(len(hist_ref.root_hists)))
    ymax = max([ hist_ref.root_flat_hist.GetBinContent(ibin + 1) + hist_ref.root_flat_hist.GetBinError(ibin + 1) for ibin in range(hist_ref.root_flat_hist.GetNbinsX()) ])
    ymin = min([ hist_ref.root_flat_hist.GetBinContent(ibin + 1) for ibin in range(hist_ref.root_flat_hist.GetNbinsX()) ])
    for hist_compare in list_hist_compare:
        ymax = max(ymax, max([ hist_compare.root_flat_hist.GetBinContent(ibin + 1) + hist_compare.root_flat_hist.GetBinError(ibin + 1) for ibin in range(hist_compare.root_flat_hist.GetNbinsX()) ]))
        ymin = min(ymin, min([ hist_compare.root_flat_hist.GetBinContent(ibin + 1) for ibin in range(hist_compare.root_flat_hist.GetNbinsX()) ]))

    Y_up = ymax * 200 if is_logY else ymax * 1.8
    hist_ref.root_flat_hist.SetAxisRange(ymin * 0.1 + 1e-06 if is_logY else 0, Y_up, 'Y')
    hist_ref.root_flat_hist.SetMarkerSize(0.5)
    hist_ref.root_flat_hist.Draw('pe')
    xAxis.SetLabelSize(0)
    xAxis.SetAxisColor(0)
    pad1.Update()
    newxAxis = [ rt.TGaxis(hist_ref.flat_xAxisLabel_low[ifun], ymin * 0.1 + 1e-06 if is_logY else 0, hist_ref.flat_xAxisLabel_up[ifun], ymin * 0.1 + 1e-06 if is_logY else 0, ('labels{}').format(ifun)) for ifun, fun in enumerate(hist_ref.flat_xAxisLabel) ]
    newxAxis_up = [ rt.TGaxis(hist_ref.flat_xAxisLabel_low[ifun], Y_up, hist_ref.flat_xAxisLabel_up[ifun], Y_up, ('labels{}').format(ifun)) for ifun, fun in enumerate(hist_ref.flat_xAxisLabel) ]
    yAxis_division = [ rt.TLine(hist_ref.flat_xAxisLabel_up[i], ymin * 0.1 + 1e-06 if is_logY else 0, hist_ref.flat_xAxisLabel_up[i], Y_up) for i in range(len(hist_ref.flat_xAxisLabel) - 1) ]
    if do_ratio:
        for axis in newxAxis:
            axis.SetLabelSize(0)

    else:
        for axis in newxAxis:
            axis.SetLabelSize(FLS)

    for iaxis, axis in enumerate(newxAxis):
        axis.SetNdivisions(205)
        if len(newxAxis) > 1 and iaxis != len(newxAxis) - 1:
            axis.ChangeLabel(-1, -1, -1, -1, -1, -1, ' ')
        axis.Draw('same')
    for axis in newxAxis_up:
        axis.SetNdivisions(205)
        axis.SetOption('-')
        axis.SetLabelSize(0)
        axis.Draw('same')

    pad1.Update()
    list_compare_legendstyle = []
    if 'Reco' in hist_ref.root_hists_name[0] and 'MC' in hist_ref.root_hists_name[0]:
        hist_ref_legendstyle = DrawPlotStyle([hist_ref.root_flat_hist], style_ref, color_ref, pad1)
    for index, hist_compare in enumerate(list_hist_compare):
        compare_legendstyle = DrawPlotStyle([hist_compare.root_flat_hist], list_style_compare[index], list_color_compare[index], pad1)
        list_compare_legendstyle.append(compare_legendstyle)
    if not ('Reco' in hist_ref.root_hists_name[0] and 'MC' in hist_ref.root_hists_name[0]):
        hist_ref_legendstyle = DrawPlotStyle([hist_ref.root_flat_hist], style_ref, color_ref, pad1)
    if len(list_hist_compare) + (hist_ref_stat != 0) < 3:
        Legend = rt.TLegend(0.9 - 0.35 / len(list_hist_compare), 0.75, 0.95, 0.93)
    else:
        Legend = rt.TLegend(0.65 - 0.35 / len(list_hist_compare), 0.75, 0.95, 0.93)
        Legend.SetNColumns(2)
    Legend.SetBorderSize(0)
    Legend.SetFillColor(0)
    Legend.SetFillStyle(0)
    Legend.SetTextFont(42)
    if do_ratio:
        Legend.SetTextSize(0.035)
    else:
        Legend.SetTextSize(0.02)
    if hist_ref_legendstyle is not None:
       Legend.AddEntry(hist_ref.root_flat_hist, legend_ref, hist_ref_legendstyle)
    for index, hist_compare in enumerate(list_hist_compare):
        if list_compare_legendstyle[index] is not None:
            Legend.AddEntry(hist_compare.root_flat_hist, list_legend_compare[index], list_compare_legendstyle[index])
    pad1.Update()
    box = create_paves(1.4817568788812e-08, 'DataPAS', CMSposX=0.17, CMSposY=0.85, prelimPosX=0.17, prelimPosY=0.8, lumiPosX=0.975, lumiPosY=0.91, alignRight=False)
    box['CMS'].Draw('same')
    box['lumi'].Draw('same')
    box['label'].Draw()
    TextPaves = []
    Interval = 0.8200000000000001 / len(text_list)
    TextCenters = np.array(range(len(text_list))) * Interval + Interval / 2 + 0.16
    if len(text_list) > 7:
        TextCenters -= 0.035
    for iText, Text in enumerate(text_list):
        TextPaves.append(rt.TPaveText(TextCenters[iText] - Interval / 2, 0.7, TextCenters[iText] + Interval / 2, 0.73, 'NDC'))
        TextPaves[iText].SetFillStyle(0)
        TextPaves[iText].SetBorderSize(0)
        TextPaves[iText].SetFillColor(0)
        TextPaves[iText].SetTextFont(42)
        if do_ratio:
            TextPaves[iText].SetTextSize(0.04)
        else:
            TextPaves[iText].SetTextSize(0.02)
        TextPaves[iText].SetTextAlign(12)
        TextPaves[iText].AddText(text_list[iText])
        TextPaves[iText].Draw('same')

    for ydiv in yAxis_division:
        ydiv.SetLineColor(1)
        ydiv.SetLineStyle(2)
        ydiv.SetLineWidth(1)
        ydiv.Draw('same')

    if do_ratio:
        rt.gStyle.SetHatchesLineWidth(1)
        c.cd()
        p1r = rt.TPad('p4', '', 0, 0, 1, 0.5)
        p1r.SetRightMargin(0.03)
        p1r.SetTopMargin(P2TM)
        p1r.SetBottomMargin(P2BM)
        p1r.SetTicks()
        p1r.Draw()
        p1r.cd()
        xmin = float(hist_ref.root_flat_hist.GetXaxis().GetXmin())
        xmax = float(hist_ref.root_flat_hist.GetXaxis().GetXmax())
        one = rt.TF1('one', '1', xmin, xmax)
        one.SetLineColor(1)
        one.SetLineStyle(2)
        one.SetLineWidth(1)
        hist_ref_ratio = hist_ref.root_flat_hist.Clone()
        if hist_ref_stat != 0:
            hist_ref_stat_ratio = hist_ref_stat.root_flat_hist.Clone()
            hist_ref_total_ratio = hist_ref_stat.root_flat_hist.Clone()
            hist_ref_stat_ratio.SetFillStyle(1001)
            hist_ref_stat_ratio.SetFillColorAlpha(rt.kGray + 2, 0.9)
            hist_ref_stat_ratio.SetLineWidth(0)
            hist_ref_stat_ratio.SetMarkerSize(0)
            hist_ref_total_ratio.SetMarkerSize(0)
            hist_ref_total_ratio.SetLineWidth(0)
            hist_ref_total_ratio.SetFillColorAlpha(rt.kGray + 1, 0.9)
            hist_ref_ratio.SetLineWidth(0)
            hist_ref_ratio.SetMarkerSize(0)
        hist_ref_ratio.SetMarkerSize(0)
        hist_ref_ratio.SetFillStyle(1001)
        hist_ref_ratio.SetFillColorAlpha(rt.kGray + 1, 0.9)
        hist_ref_ratio.SetLineWidth(0)
        list_hist_compare_ratio = []
        ratio_min = 1.0
        ratio_max = 1.0
        for hist_compare in list_hist_compare:
            list_hist_compare_ratio.append(hist_compare.root_flat_hist.Clone())

        for ibin in range(hist_ref.root_flat_hist.GetNbinsX()):
            ndata = hist_ref.root_flat_hist.GetBinContent(ibin + 1)
            edata = hist_ref.root_flat_hist.GetBinError(ibin + 1)
            hist_ref_ratio.SetBinContent(ibin + 1, 1.0)
            hist_ref_ratio.SetBinError(ibin + 1, edata / ndata if ndata > 0 else 0)
            for index, hist_compare_ratio in enumerate(list_hist_compare_ratio):
                hist_compare_ratio.SetBinContent(ibin + 1, hist_compare_ratio.GetBinContent(ibin + 1) / ndata if ndata > 0 else 0)
                hist_compare_ratio.SetBinError(ibin + 1, hist_compare_ratio.GetBinError(ibin + 1) / ndata if ndata > 0 else 0)

            if hist_ref_stat != 0:
                hist_ref_total_ratio.SetBinContent(ibin + 1, 1.0)
                hist_ref_total_ratio.SetBinError(ibin + 1, edata / ndata if ndata > 0 else 0)
                hist_ref_stat_ratio.SetBinContent(ibin + 1, hist_ref_stat.root_flat_hist.GetBinContent(ibin + 1) / ndata if ndata > 0 else 0)
                hist_ref_stat_ratio.SetBinError(ibin + 1, hist_ref_stat.root_flat_hist.GetBinError(ibin + 1) / ndata if ndata > 0 else 0)

        hist_ref_ratio.SetAxisRange(max(0, ratio_min - 0.3), min(2, ratio_max + 0.3), 'Y')
        hist_ref_ratio.SetTitle('')
        hist_ref_ratio.GetXaxis().SetTitle(title)
        hist_ref_ratio.GetXaxis().SetTitleSize(RTSX)
        hist_ref_ratio.GetXaxis().SetTitleOffset(RTOX)
        hist_ref_ratio.GetXaxis().SetLabelSize(0)
        hist_ref_ratio.GetXaxis().SetAxisColor(0)
        hist_ref_ratio.GetXaxis().SetLabelOffset(0.02)
        hist_ref_ratio.GetYaxis().SetTitleSize(RTSY)
        hist_ref_ratio.GetYaxis().SetLabelSize(RLSY)
        hist_ref_ratio.GetYaxis().SetTitleOffset(0.2 + 0.15 / len(hist_ref.root_hists))
        hist_ref_ratio.GetYaxis().SetTitle('            ' + label_ratio)
        if len(label_ratio) > 10:
            hist_ref_ratio.GetYaxis().SetTitleSize(RTSY * 0.6)
            hist_ref_ratio.GetYaxis().SetTitleOffset(0.5 + 0.22 / len(hist_ref.root_hists))
        hist_ref_ratio.GetYaxis().SetDecimals(1)
        hist_ref_ratio.GetYaxis().SetMaxDigits(3)
        hist_ref_ratio.GetYaxis().SetNdivisions(4, 2, 0)
        p1r.cd()
        hist_ref_ratio.Draw('e2same')
        newxAxisratio = [ rt.TGaxis(hist_ref.flat_xAxisLabel_low[ifun], max(0, ratio_min - 0.3), hist_ref.flat_xAxisLabel_up[ifun], max(0, ratio_min - 0.3), ('labels{}').format(ifun)) for ifun, fun in enumerate(hist_ref.flat_xAxisLabel) ]
        for iaxis, axis in enumerate(newxAxisratio):
            axis.SetLabelSize(0.05)
            axis.SetLabelOffset(0.02)
            axis.SetNdivisions(205)
            if len(newxAxisratio) > 1 and iaxis != len(newxAxisratio) - 1:
                axis.ChangeLabel(-1, -1, -1, -1, -1, -1, ' ')
            axis.Draw('same')

        newxAxisratio_up = [ rt.TGaxis(hist_ref.flat_xAxisLabel_low[ifun], min(2, ratio_max + 0.3), hist_ref.flat_xAxisLabel_up[ifun], min(2, ratio_max + 0.3), ('labels{}').format(ifun)) for ifun, fun in enumerate(hist_ref.flat_xAxisLabel) ]
        for iaxis, axis in enumerate(newxAxisratio_up):
            axis.SetNdivisions(205)
            axis.SetOption('-')
            axis.SetLabelSize(0)
            axis.Draw('same')

        yAxisratio_division = [ rt.TLine(hist_ref.flat_xAxisLabel_up[i], max(0, ratio_min - 0.3), hist_ref.flat_xAxisLabel_up[i], min(2, ratio_max + 0.3)) for i in range(len(hist_ref.flat_xAxisLabel) - 1) ]
        for ydivratio in yAxisratio_division:
            ydivratio.SetLineColor(1)
            ydivratio.SetLineStyle(2)
            ydivratio.SetLineWidth(1)
            ydivratio.Draw('same')

        p1r.Update()
        if hist_ref_stat != 0:
            p1r.cd()
            hist_ref_total_ratio.Draw('e2same')
            hist_ref_stat_ratio.Draw('e2same')
            p1r.Update()
            pad1.cd()
            Legend.AddEntry(hist_ref_stat_ratio, 'Stat. unc', 'f')
            Legend.AddEntry(hist_ref_total_ratio, 'Sys.+Stat. unc', 'f')
            p1r.cd()
        if 'statonly' in hist_ref.root_hists[0].GetName():
            pad1.cd()
            Legend.AddEntry(hist_ref_ratio, 'Stat. unc', 'f')
            p1r.cd()
        for index, hist_compare_ratio in enumerate(list_hist_compare_ratio):
            if list_style_compare[index] == 'filled':
                style_ratio = 'fillederror'
            else:
                style_ratio = list_style_compare[index]
            DrawPlotStyle([hist_compare_ratio], style_ratio, list_color_compare[index], p1r)

        one.Draw('same')
        p1r.Update()
    pad1.cd()
    Legend.Draw()
    pad1.Update()
    c.Update()
    c.SaveAs(output_path + '.pdf')
    c.SaveAs(output_path + '.png')
    return
