import numpy as np
import json
from argparse import ArgumentParser
import os,ast
import time
import math
import sys,stat
import ROOT as rt
import itertools
import  tdrstyle
import CMS_lumi
from Plotting_cfg import *
import h5py
import pandas as pd
import itertools
from math import sqrt
from array import array
import glob
from sklearn import metrics
from scipy.stats import chi2 as CHI2
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from plot_utils import *
from unfold_utils import *
def GOF(HistList1,HistList2):
  assert len(HistList1)==len(HistList2)
  wchi2=0.
  ndof=0.
  wchi2_hist2unc=0.
  for (hist1,hist2) in zip(HistList1,HistList2):
    assert hist1.GetNbinsX()==hist2.GetNbinsX()
    wchi2 += sum([(hist1.GetBinContent(ibin+1)-hist2.GetBinContent(ibin+1))**2/(hist1.GetBinError(ibin+1)**2+hist2.GetBinError(ibin+1)**2) for ibin in range(hist1.GetNbinsX())])
    wchi2_hist2unc += sum([(hist1.GetBinContent(ibin+1)-hist2.GetBinContent(ibin+1))**2/(hist2.GetBinError(ibin+1)**2) for ibin in range(hist1.GetNbinsX())])
    ndof += hist1.GetNbinsX()
  ndof -= 1
  p= CHI2.sf(wchi2,ndof)
  p_hist2unc= CHI2.sf(wchi2_hist2unc,ndof)
  return wchi2,wchi2_hist2unc,ndof,wchi2/ndof,wchi2_hist2unc/ndof,p,p_hist2unc

def GetRefoldIter(File):
  names=[key.GetName() for key in File.GetListOfKeys() if ("HistRecoInclusive_MC" in key.GetName() and "iter" in key.GetName())]
  #extract heads and remove duplicates
  heads=list(dict.fromkeys([name.split("iter")[0] for name in names]))
  print heads
  heads_sort=sorted(heads,key=lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])
  print heads_sort
  if any("iter0" in s for s in names):
    start_iter=0
  else:
    start_iter=1
  names_rank=[[head+"iter"+str(i) for head in heads_sort] for i in range(start_iter,len(names)/len(heads_sort)+start_iter)]
  print names_rank
  names_data=[key.GetName() for key in File.GetListOfKeys() if "HistRecoInclusive_Pseudodata" in key.GetName()]
  if len(names_data)>0:
    names_data_sort=sorted(names_data,key=lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])
  else:
    names_data=[key.GetName() for key in File.GetListOfKeys() if "HistRecoInclusive_Data" in key.GetName()]
    names_data_sort=sorted(names_data,key=lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])
  names_MC=[key.GetName() for key in File.GetListOfKeys() if "HistRecoInclusive_MC_" in key.GetName() and "iter" not in key.GetName()]
  names_MC_sort=sorted(names_MC,key=lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])
  return names_rank,names_data_sort,names_MC_sort

def GetUnfoldIter(File):
  names=[key.GetName() for key in File.GetListOfKeys() if ("HistGenInclusive_MC" in key.GetName() and "iter" in key.GetName())]
  heads=list(dict.fromkeys([name.split("iter")[0] for name in names]))
  heads_sort=sorted(heads,key=lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])
  if any("iter0" in s for s in names):
    start_iter=0
  else:
    start_iter=1
  names_rank=[[head+"iter"+str(i) for head in heads_sort] for i in range(start_iter,len(names)/len(heads_sort)+start_iter)]
  names_MC=[key.GetName() for key in File.GetListOfKeys() if "HistGenInclusive_MC" in key.GetName() and "iter" not in key.GetName()]
  names_MC_sort=sorted(names_MC,key=lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])
  return names_rank,names_MC_sort

def GetPseudoDataTruth(File):
  names=[key.GetName() for key in File.GetListOfKeys() if "HistGenInclusive_Pseudodata" in key.GetName()]
  names_sort=None
  if len(names)>0:
    names_sort=sorted(names,key=lambda s: np.array(re.findall(r'\d+', s),dtype=int)[-1])
  return names_sort


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('--input',default='results_finebin_v7_MCCP1ES_CP5sys_trksys_1d_optimize/unfold_nparticle_eta2p4pt05_pur_1d_spherocity_eta2p4pt05_pur_1d_nominal_optimize_omnifold.root',help='The input root file containing the results in iterations')
    parser.add_argument('--output',default='results_finebin_v7_MCCP1ES_CP5sys_trksys_1d_optimize/iter_nparticle_eta2p4pt05_pur_1d_spherocity_eta2p4pt05_pur_1d_nominal_omnifold.root',help='The output root file containing the refolfing chi^2 w.r.t. iterations')
    parser.add_argument('--config',default="Config_plot/Config_sph_1d_v7_MCCP1ES_CP5sys.json",help="The configration file including the unfolding setup")
    parser.add_argument('--plot',action="store_true",default=False)
    parser.add_argument('--plotdir',default="results_finebin_v7_MCCP1ES_CP5sys_trksys_1d_optimize/plots_optimize")
    parser.add_argument('--sysreweight',action="store_true",default=False)
    args = parser.parse_args()
    with open(args.config, 'r') as configjson:
        config = json.load(configjson)
    with open(config["varunfold"], 'r') as fjson:
        info_var = json.load(fjson)
    FineBin = "binedgesreco" in info_var[config["var1"]].keys() and "binedgesreco" in info_var[config["var2"]].keys()
    TextListReco=[]
    TextListGen=[]
    if FineBin:
      TextListReco = [str(info_var[config["var1"]]["binedgesreco"][i])+"#leq"+info_var[config["var1"]]["reco_shortname"]+"<"+str(info_var[config["var1"]]["binedgesreco"][i+1]) for i in range(info_var[config["var1"]]["nbinsreco"])]
      TextListGen = [str(info_var[config["var1"]]["binedgesgen"][i])+"#leq"+info_var[config["var1"]]["gen_shortname"]+"<"+str(info_var[config["var1"]]["binedgesgen"][i+1]) for i in range(info_var[config["var1"]]["nbinsgen"])]+(["background"] if config["addbkg"] else [])
    else:
      TextListReco = [str(info_var[config["var1"]]["minreco"]+(info_var[config["var1"]]["maxreco"]-info_var[config["var1"]]["minreco"])/info_var[config["var1"]]["nbinsreco"]*i)+"#leq"+info_var[config["var1"]]["reco_shortname"]+"<"+str(info_var[config["var1"]]["minreco"]+(info_var[config["var1"]]["maxreco"]-info_var[config["var1"]]["minreco"])/info_var[config["var1"]]["nbinsreco"]*(i+1)) for i in range(info_var[config["var1"]]["nbinsreco"])]
      TextListGen = [str(info_var[config["var1"]]["mingen"]+(info_var[config["var1"]]["maxgen"]-info_var[config["var1"]]["mingen"])/info_var[config["var1"]]["nbinsgen"]*i)+"#leq"+info_var[config["var1"]]["gen_shortname"]+"<"+str(info_var[config["var1"]]["mingen"]+(info_var[config["var1"]]["maxgen"]-info_var[config["var1"]]["mingen"])/info_var[config["var1"]]["nbinsgen"]*(i+1)) for i in range(info_var[config["var1"]]["nbinsgen"])]+(["background"] if config["addbkg"] else [])
    f=rt.TFile(args.input,"READ")

    names_refold,name_data,name_MCreco=GetRefoldIter(f)
    names_unfold,name_MC=GetUnfoldIter(f)
    names_psedodata_truth=GetPseudoDataTruth(f)
    fout=rt.TFile(args.output,"RECREATE")
    hist_chi2_iter_dataMCunc=rt.TH1F("Refold_chi2_dataMCunc","Refold #chi^{2}/n.d.o.f",len(names_refold),0,len(names_refold))
    hist_chi2_iter_dataMCunc.GetXaxis().SetTitle("iterations")
    hist_chi2_iter_dataMCunc.GetYaxis().SetTitle("#chi^{2}/n.d.o.f")
    hist_chi2_iter_dataunc=rt.TH1F("Refold_chi2_dataunc","Refold #chi^{2}/n.d.o.f",len(names_refold),0,len(names_refold))
    hist_chi2_iter_dataunc.GetXaxis().SetTitle("iterations")
    hist_chi2_iter_dataunc.GetYaxis().SetTitle("#chi^{2}/n.d.o.f")

    hist_chi2_iter_MC=rt.TH1F("Unfold_chi2_MC","Unfold #chi^{2}/n.d.o.f",len(names_unfold),0,len(names_unfold))
    hist_chi2_iter_MC.GetXaxis().SetTitle("iterations")
    hist_chi2_iter_MC.GetYaxis().SetTitle("#chi^{2}/n.d.o.f")

    hist_chi2_iter_MC_MCunfoldunc=rt.TH1F("Unfold_chi2_MC_MCunfoldunc","Unfold #chi^{2}/n.d.o.f",len(names_unfold),0,len(names_unfold))
    hist_chi2_iter_MC_MCunfoldunc.GetXaxis().SetTitle("iterations")
    hist_chi2_iter_MC_MCunfoldunc.GetYaxis().SetTitle("#chi^{2}/n.d.o.f")

    hist_p_iter_MC_MCunfoldunc=rt.TH1F("Unfold_p_MC_MCunfoldunc","Unfold p",len(names_unfold),0,len(names_unfold))
    hist_p_iter_MC_MCunfoldunc.GetXaxis().SetTitle("iterations")
    hist_p_iter_MC_MCunfoldunc.GetYaxis().SetTitle("p-value")

    if not(names_psedodata_truth is None):
      hist_chi2_iter_truth_dataunc=rt.TH1F("Unfold_chi2_pseudodatatruth_dataunc","Unfold #chi^{2}/n.d.o.f",len(names_unfold),0,len(names_unfold))
      hist_chi2_iter_truth_dataunc.GetXaxis().SetTitle("iterations")
      hist_chi2_iter_truth_dataunc.GetYaxis().SetTitle("#chi^{2}/n.d.o.f")
      hist_chi2_iter_truth_dataMCunc=rt.TH1F("Unfold_chi2_pseudodatatruth_dataMCunc","Unfold #chi^{2}/n.d.o.f",len(names_unfold),0,len(names_unfold))
      hist_chi2_iter_truth_dataMCunc.GetXaxis().SetTitle("iterations")
      hist_chi2_iter_truth_dataMCunc.GetYaxis().SetTitle("#chi^{2}/n.d.o.f")


    _,_,_,wchi2_per_ndof_MCdatareco,wchi2_hist2unc_per_ndof_MCdatareco,p_MCdatareco,_ = GOF([f.Get(name) for name in name_MCreco],[f.Get(name) for name in name_data])
    if not(names_psedodata_truth is None):
      _,_,_,wchi2_per_ndof_MCdatagen,wchi2_hist2unc_per_ndof_MCdatagen,_,_ = GOF([f.Get(name) for name in name_MC],[f.Get(name) for name in names_psedodata_truth])

    hist_list_data=hist_list("Data")
    hist_list_data.root_hists_name=name_data
    hist_list_data.read_hist_from_file(f)
    hist_list_data.divide_by_bin_width()
    hist_list_data.flatten_hist()
    hist_list_MCgeninclusive=hist_list("MCGenInclusive")
    hist_list_MCgeninclusive.root_hists_name=name_MC
    hist_list_MCgeninclusive.read_hist_from_file(f)
    hist_list_MCgeninclusive.divide_by_bin_width()
    hist_list_MCgeninclusive.flatten_hist()
    hist_list_MCrecoinclusive=hist_list("MCRecoInclusive")
    hist_list_MCrecoinclusive.root_hists_name=name_MCreco
    hist_list_MCrecoinclusive.read_hist_from_file(f)
    hist_list_MCrecoinclusive.divide_by_bin_width()
    hist_list_MCrecoinclusive.flatten_hist()
    if not(names_psedodata_truth is None):
      hist_list_pseudodatatruthinclusive=hist_list("PseudodataTruthInclusive")
      hist_list_pseudodatatruthinclusive.root_hists_name=names_psedodata_truth
      hist_list_pseudodatatruthinclusive.read_hist_from_file(f)
      hist_list_pseudodatatruthinclusive.divide_by_bin_width()
      hist_list_pseudodatatruthinclusive.flatten_hist()
    else:
      hist_list_pseudodatatruthinclusive=None

    Config={}
    if args.plot:
      Config["Data"]={
        "hist":hist_list_data,
        "stat":0,
        "color":1,
        "style":"cross",
        "legend":("Data" if names_psedodata_truth is None else "Pseudo-data") if not args.sysreweight else "sys variation: "+config["syslegend"][0]
      }
      Config["MCGenInclusive"]={
        "hist":hist_list_MCgeninclusive,
        "stat":0,
        "color":rt.kAzure-2,
        "style":"fillederror",
        "legend":str(config["MClegend"])
      }
      Config["MCRecoInclusive"]={
        "hist":hist_list_MCrecoinclusive,
        "stat":0,
        "color":rt.kAzure-2,
        "style":"fillederror",
        "legend":str(config["MClegend"])
      }
      Config["PseudodataTruthInclusive"]={
        "hist":hist_list_pseudodatatruthinclusive,
        "stat":0,
        "color":800,
        "style":"triangle",
        "legend":"Pseudo-data truth" if not args.sysreweight else "sys variation: "+config["syslegend"][0]
      }
    PlotLists={}
    def PlotConfig(PlotName):
      comparehists=[c["hist"] for c in PlotName["compare"] if c["hist"] is not None]
      if len(comparehists)==0:
        return
      comparecolors=[c["color"] for c in PlotName["compare"] if c["hist"] is not None]
      comparestyles=[c["style"] for c in PlotName["compare"] if c["hist"] is not None]
      comparelegends=[c["legend"] for c in PlotName["compare"] if c["hist"] is not None]
      path=str(args.plotdir+"/"+PlotName["name"]+"_"+config["var1"]+"_"+config["var2"])
      os.system("mkdir -p "+args.plotdir)
      if "reco" in PlotName["name"] or "refold" in PlotName["name"]:
        axis_title=info_var[config["var2"]]["reco_name"]
      else:
        axis_title=info_var[config["var2"]]["gen_name"]
      plot_flat_hists(PlotName["ref"]["hist"], comparehists, PlotName["ref"]["legend"], comparelegends, title=axis_title, is_logY=1, do_ratio=1, output_path=path+"_logy", hist_ref_stat=PlotName["ref"]["stat"], text_list=(TextListReco if ("reco" in PlotName["name"] or "refold" in PlotName["name"]) else TextListGen), style_ref=PlotName["ref"]["style"], color_ref=PlotName["ref"]["color"], list_style_compare=comparestyles, list_color_compare=comparecolors, labelY='Normalized Events/Bin Width', label_ratio=PlotName["ratio"])
      return path+"_logy.png"

    print(name_data)
    for i,name_refold in enumerate(names_refold):
      print(name_refold)
      _,_,_,wchi2_per_ndof,wchi2_hist2unc_per_ndof,_,_ = GOF([f.Get(name) for name in name_refold],[f.Get(name) for name in name_data])
      hist_chi2_iter_dataMCunc.SetBinContent(i+1,wchi2_per_ndof)
      hist_chi2_iter_dataMCunc.SetBinError(i+1,0)
      hist_chi2_iter_dataunc.SetBinContent(i+1,wchi2_hist2unc_per_ndof)
      hist_chi2_iter_dataunc.SetBinError(i+1,0)

    for i,name_unfold in enumerate(names_unfold):
      _,_,_,wchi2_per_ndof,wchi2_hist2unc_per_ndof,p,_ = GOF([f.Get(name) for name in name_unfold],[f.Get(name) for name in name_MC])
      iter_index=re.search("iter(\d+)",name_unfold[0]).group(1)
      hist_chi2_iter_MC.SetBinContent(i+1,wchi2_hist2unc_per_ndof)
      hist_chi2_iter_MC.SetBinError(i+1,0)
      hist_chi2_iter_MC_MCunfoldunc.SetBinContent(i+1,wchi2_per_ndof)
      hist_chi2_iter_MC_MCunfoldunc.SetBinError(i+1,0)
      hist_p_iter_MC_MCunfoldunc.SetBinContent(i+1,p)
      hist_p_iter_MC_MCunfoldunc.SetBinError(i+1,0)
      if not(names_psedodata_truth is None):
        _,_,_,wchi2_per_ndof,wchi2_hist2unc_per_ndof,_,_ = GOF([f.Get(name) for name in name_unfold],[f.Get(name) for name in names_psedodata_truth])
        hist_chi2_iter_truth_dataunc.SetBinContent(i+1,wchi2_hist2unc_per_ndof)
        hist_chi2_iter_truth_dataunc.SetBinError(i+1,0)
        hist_chi2_iter_truth_dataMCunc.SetBinContent(i+1,wchi2_per_ndof)
        hist_chi2_iter_truth_dataMCunc.SetBinError(i+1,0)
      if args.plot:
        color=4
        legend="MLE refold" if not args.sysreweight else "MLE reweight"
        tag="MLE"
        color_unfold=4
        legend_unfold="MLE unfold" if not args.sysreweight else "MLE reweight"
        if "omnifold" in args.input:
          color=608
          legend="omnifold refold" if not args.sysreweight else "omnifold reweight"
          tag="omnifold"
          color_unfold=6
          legend_unfold="omnifold unfold" if not args.sysreweight else "omnifold reweight"
        elif "multifold" in args.input:
          color=812
          legend="multifold refold" if not args.sysreweight else "multifold reweight"
          tag="multifold"
          color_unfold=8
          legend_unfold="multifold unfold" if not args.sysreweight else "multifold reweight"
        elif "unifold" in args.input:
          color=46
          legend="unifold refold" if not args.sysreweight else "unifold reweight"
          tag="unifold"
          color_unfold=2
          legend_unfold="unifold unfold" if not args.sysreweight else "unifold reweight"
        hist_list_refold=hist_list("Refold_iter"+iter_index)
        hist_list_refold.root_hists_name=names_refold[i]
        hist_list_refold.read_hist_from_file(f)
        hist_list_refold.divide_by_bin_width()
        hist_list_refold.flatten_hist() 
        hist_list_unfold=hist_list("Unfold_iter"+iter_index)
        hist_list_unfold.root_hists_name=names_unfold[i]
        hist_list_unfold.read_hist_from_file(f)
        hist_list_unfold.divide_by_bin_width()
        hist_list_unfold.flatten_hist() 
        Config["Refold_iter"+iter_index]={
          "hist":hist_list_refold,
          "stat":0,
          "color":color,
          "style":"filled",
          "legend":legend
        }
        Config["Unfold_iter"+iter_index]={
          "hist":hist_list_unfold,
          "stat":0,
          "color":color_unfold,
          "style":"marker",
          "legend":legend_unfold
        }
        PlotLists["Refoldcompare_iter"+iter_index]={
          "ref":Config["Data"],
          "compare":[Config["Refold_iter"+iter_index],Config["MCRecoInclusive"]],
          "name":"data_refold_"+tag+"_iter"+iter_index,
          "ratio":"Refold / data" if not args.sysreweight else "Reweight / sys. var."
        }
        PlotLists["Unfoldcompare_iter"+iter_index]={
          "ref":Config["MCGenInclusive"],
          "compare":[Config["Unfold_iter"+iter_index]],
          "name":"MC_unfold_"+tag+"_iter"+iter_index,
          "ratio":"Unfold / MC" if not args.sysreweight else "Reweight / MC"
        }
        PlotLists["Unfoldcomparepseudodata_iter"+iter_index]={
          "ref":Config["PseudodataTruthInclusive"],
          "compare":[Config["Unfold_iter"+iter_index],Config["MCGenInclusive"]],
          "name":"pseudodatatruth_unfold_"+tag+"_iter"+iter_index,
          "ratio":"Unfold / Truth" if not args.sysreweight else "Reweight / sys. var."
        }
        PlotConfig(PlotLists["Refoldcompare_iter"+iter_index])
        print("plotting")
        PlotConfig(PlotLists["Unfoldcompare_iter"+iter_index])
        if not(names_psedodata_truth is None):
          PlotConfig(PlotLists["Unfoldcomparepseudodata_iter"+iter_index])
    hist_chi2_iter_dataMCunc.Write()
    hist_chi2_iter_dataunc.Write()
    hist_chi2_iter_MC.Write()
    hist_chi2_iter_MC_MCunfoldunc.Write()
    hist_p_iter_MC_MCunfoldunc.Write()
    line_reco_dataMCunc=rt.TLine(0,wchi2_per_ndof_MCdatareco,len(names_refold),wchi2_per_ndof_MCdatareco)
    line_reco_dataMCunc.SetLineColor(rt.kRed)
    line_reco_dataMCunc.Write("Chi2dataMC_reco_dataMCunc")
    line_reco_dataunc=rt.TLine(0,wchi2_hist2unc_per_ndof_MCdatareco,len(names_refold),wchi2_hist2unc_per_ndof_MCdatareco)
    line_reco_dataunc.SetLineColor(rt.kRed)
    line_reco_dataunc.Write("Chi2dataMC_reco_dataunc")

    line_reco_dataMCunc_p=rt.TLine(0,p_MCdatareco,len(names_refold),p_MCdatareco)
    line_reco_dataMCunc_p.SetLineColor(rt.kRed)
    line_reco_dataMCunc_p.Write("p_dataMC_reco_dataMCunc")

    if not(names_psedodata_truth is None):
      hist_chi2_iter_truth_dataunc.Write()
      hist_chi2_iter_truth_dataMCunc.Write()
      line_gen_dataunc=rt.TLine(0,wchi2_hist2unc_per_ndof_MCdatagen,len(names_unfold),wchi2_hist2unc_per_ndof_MCdatagen)
      line_gen_dataunc.SetLineColor(rt.kRed)
      line_gen_dataunc.Write("Chi2pseudodataMC_gen_pseudodataunc")
      line_gen_dataMCunc=rt.TLine(0,wchi2_per_ndof_MCdatagen,len(names_unfold),wchi2_per_ndof_MCdatagen)
      line_gen_dataMCunc.SetLineColor(rt.kRed)
      line_gen_dataMCunc.Write("Chi2pseudodataMC_gen_pseudodataMCunc")

    fout.Close()

