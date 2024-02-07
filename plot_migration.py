import numpy as np
import json
from argparse import ArgumentParser
import os,ast
import time
import math
import sys,stat
import ROOT
import itertools
import h5py
import pandas as pd
from array import array
import math

ROOT.gROOT.SetBatch(ROOT.kTRUE)
if __name__=="__main__":

    parser = ArgumentParser()
    parser.add_argument('--inputref', default="results_finebin_v7_MCA3P_sysweightCH3genweightA3P_1d_ensemble4_v3/unfold_nparticle_spherocity_nominal_optimize_multifold_step1_debug.root", help="The input root file for plotting (histref)")
    parser.add_argument('--inputcompare', default="results_finebin_v7_MCA3P_sysweightCH3genweightA3P_1d_ensemble4_v3/unfold_nparticle_spherocity_nominal_optimize_multifold_step1_debug.root", help="The input root file for plotting (histcompare)")
    parser.add_argument('--output', default = "results_finebin_v7_MCA3P_sysweightCH3genweightA3P_1d_ensemble4_v3/plots_mig/mig_compare",help="output name")
    parser.add_argument('--histref',default="HistMig_MC_multifold_gen_nch3.0-140.0_reco_ntrk3.0-160.0_iter0",help="The name of the histogram for plotting")
    parser.add_argument('--histcompare',default=None,help="If not provided, then plot the histogram histref, otherwise plot the relative difference between the histref and histcompare")
    parser.add_argument('--title',default="Migration matrix of MC (spherocity)",help="The title of the plot")
    parser.add_argument('--row',default=1,type=int,help="number of rows (#reco coarse bins)")
    parser.add_argument('--column',default=1,type=int,help="number of column (#gen coarse bins)")
    parser.add_argument('--significance',action="store_true",default=False)
    parser.add_argument('--range',default=None,type=int,help="range of the z axis")
    args = parser.parse_args()
    #ROOT.gInterpreter.Declare("TH2F* recast(const char *path, const char *funcname) { return (TH2F *) TFile(path).Get(funcname)->Clone();}")
    fref = ROOT.TFile(args.inputref,"read")
    histref_name_list = args.histref.split(",")
    histref_list = []
    for histref_name in histref_name_list:
      histref = fref.Get(histref_name)
      #print(histref)
      histref_list.append(histref) 
    #print(histref_list)
    os.makedirs(os.path.dirname(args.output),exist_ok=True)
    if args.histcompare:
      fcompare = ROOT.TFile(args.inputcompare,"read")
      histcompare_name_list = args.histcompare.split(",")
      histcompare_list = []
      for histcompare_name in histcompare_name_list:
        histcompare = fcompare.Get(histcompare_name)
        #print(histcompare)
        histcompare_list.append(histcompare)
      #print(histcompare_list)
      histdiff_list = []
      histsignificance_list = []
      for i,histcompare in enumerate(histcompare_list):
        histdiff_list.append(histcompare.Clone())
        histdiff_list[i].Add(histcompare,histref_list[i],1,-1)
        histsignificance_list.append(histdiff_list[i].Clone())
        for ibin1 in range(histsignificance_list[i].GetNbinsX()):
          for ibin2 in range(histsignificance_list[i].GetNbinsY()):
            #histsignificance_list[i].SetBinContent(ibin1+1,ibin2+1,(histsignificance_list[i].GetBinContent(ibin1+1,ibin2+1)/np.sqrt(histsignificance_list[i].GetBinError(ibin1+1,ibin2+1)) if histsignificance_list[i].GetBinError(ibin1+1,ibin2+1)>0 else 0))
            histsignificance_list[i].SetBinContent(ibin1+1,ibin2+1,(histsignificance_list[i].GetBinContent(ibin1+1,ibin2+1)/histsignificance_list[i].GetBinError(ibin1+1,ibin2+1) if histsignificance_list[i].GetBinError(ibin1+1,ibin2+1)>0 else 0))
        histdiff_list[i].Divide(histref_list[i])
    if args.histcompare:
      if args.significance:
        histplot_list = histsignificance_list
      else:
        histplot_list = histdiff_list
    else:
      histplot_list = histref_list

    #if not args.histcompare:
    #  for histplot in histplot_list:
    #    for ibinx in range(histplot.GetNbinsX()):
    #      for ibiny in range(histplot.GetNbinsY()):
    #        histplot.SetBinError(ibinx+1,ibiny+1,histplot.GetBinError(ibinx+1,ibiny+1)/histplot.GetYaxis().GetBinWidth(ibiny+1))
    #        histplot.SetBinContent(ibinx+1,ibiny+1,histplot.GetBinContent(ibinx+1,ibiny+1)/histplot.GetYaxis().GetBinWidth(ibiny+1))

    if args.histcompare:
      red   = array('d',[ 0.80, 0.90, 1.00, 0.60, 0.02, ])
      green = array('d',[ 0.20, 0.80, 1.00, 0.80, 0.20, ])
      blue  = array('d',[ 0.10, 0.60, 1.00, 0.90, 0.65, ])
    else:
      red   = array('d',[ 1.00, 1., ])
      green = array('d',[ 1.00, 0., ])
      blue  = array('d',[ 1.00, 0., ])
    if args.histcompare:
      stops = array('d',[i/(len(red)-1.) for i in range(0,len(red))])
      FI = ROOT.TColor.CreateGradientColorTable(len(red), stops, red, green, blue, 100)
      kMyTemperature = array('i',[ FI+i for i in range(100)])
      ROOT.gStyle.SetPalette(100,kMyTemperature)
      ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
    ROOT.gStyle.SetOptStat(0)
    limit=0
    if args.histcompare:
      if args.significance:
        for histplot in histplot_list:
          print("maximum: ", histplot.GetMaximum(), "minimum: ", histplot.GetMinimum())
          if args.range: limit = args.range
          else:
            limit = math.ceil(max(abs(histplot.GetMaximum()),abs(histplot.GetMinimum())))
          print("setting limit to: ",limit)
          histplot.GetZaxis().SetRangeUser(-limit,limit)
      else:
        for histplot in histplot_list:
          print("maximum: ", histplot.GetMaximum(), "minimum: ", histplot.GetMinimum())
          if args.range: limit = args.range
          else:
            limit = math.ceil(max(abs(histplot.GetMaximum()),abs(histplot.GetMinimum())))
          print("setting limit to: ",limit)
          histplot.GetZaxis().SetRangeUser(-limit,limit)
    #else:
    #  for histplot in histplot_list:
    #    histplot.GetZaxis().SetRangeUser(0,1)

    title_list = args.title.split(",")
    c = ROOT.TCanvas()
    c.Divide(args.column,args.row,0,0)
    ROOT.gStyle.SetPaintTextFormat("1.3f")
    for i,histplot in enumerate(histplot_list):
      histplot.GetXaxis().SetTitle("Gen")
      histplot.GetYaxis().SetTitle("Reco")
      histplot.SetTitle(title_list[i])
      #histplot.SetMarkerSize(0.2)
      c.cd(i+1)
      if not args.histcompare: ROOT.gPad.SetLogz()
      if i%args.column == args.column-1:
        ROOT.gPad.SetRightMargin(0.1)
      histplot.Draw("colz")
      c.Update()
      #histplot.Draw("colz ERROR TEXT45")
    if not args.histcompare: args.output+="_log"
    c.SaveAs(args.output+".pdf")
    c.SaveAs(args.output+".png")
    sys.exit(limit)


