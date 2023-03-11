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
import itertools
from math import sqrt
from array import array
import glob
from unfold_utils import *

def fill_hist_lists(dataset,var1_config,var2_config,edges_gen,edges_reco,source,genWeight="",from_root=True,weight_array=None,store_mig=False,tag="", reco_only=False):
  hists = {}

  reco_inclusive=hist_list("HistRecoInclusive_"+dataset,tag)
  reco_inclusive.read_settings_from_config_dim1(var1_config,isgen=False)
  reco_inclusive.read_settings_from_config_dim2(var2_config,isgen=False)
  reco_inclusive.bin_edges_dim2 = edges_reco
  reco_inclusive.cut = cuts[CutType.PassReco]
  hists["reco_inclusive"] = reco_inclusive

  if not reco_only:
    gen_passreco=hist_list("HistGen_"+dataset,tag)
    gen_passreco.read_settings_from_config_dim1(var1_config,isgen=True)
    gen_passreco.read_settings_from_config_dim2(var2_config,isgen=True)
    gen_passreco.bin_edges_dim2 = edges_gen
    gen_passreco.cut = cuts[CutType.PassReco_PassGen]
    hists["gen_passreco"] = gen_passreco

    gen_inclusive=hist_list("HistGenInclusive_"+dataset,tag)
    gen_inclusive.read_settings_from_config_dim1(var1_config,isgen=True)
    gen_inclusive.read_settings_from_config_dim2(var2_config,isgen=True)
    gen_inclusive.bin_edges_dim2 = edges_gen
    gen_inclusive.cut = cuts[CutType.PassGen]
    hists["gen_inclusive"] = gen_inclusive

    reco_passgen=hist_list("HistReco_"+dataset,tag)
    reco_passgen.read_settings_from_config_dim1(var1_config,isgen=False)
    reco_passgen.read_settings_from_config_dim2(var2_config,isgen=False)
    reco_passgen.bin_edges_dim2 = edges_reco
    reco_passgen.cut = cuts[CutType.PassReco_PassGen]
    hists["reco_passgen"] = reco_passgen

  for _, hist in hists.items():
    hist.fill_root_hists_name()
    print "histograms:", hist.root_hists_name
    hist.fill_hist(source, from_root, weightarray=weight_array,genWeight=genWeight)
  print "filled"

  if store_mig:
    mig = [[hist_list("HistMig_"+dataset+"_"+var1_config["gen_key"]+str(var1_config["binedgesgen"][i])+"-"+str(var1_config["binedgesgen"][i+1])+"_"+var1_config["reco_key"]+str(var1_config["binedgesreco"][j])+"-"+str(var1_config["binedgesreco"][j+1])+tag) for j in range(len(var1_config["binedgesreco"])-1)] for i in range(len(var1_config["binedgesgen"])-1)]
    for i in range(len(var1_config["binedgesgen"])-1):
      for j in range(len(var1_config["binedgesreco"])-1):
        mig[i][j].read_settings_from_config_dim1(var2_config,isgen=True)
        mig[i][j].read_settings_from_config_dim2(var2_config,isgen=False)
        mig[i][j].root_cut=root_cuts[CutType.PassReco_PassGen]+"*({var1gen}>={var1gen_low})*({var1gen}<{var1gen_high})*({var1reco}>={var1reco_low})*({var1reco}<{var1reco_high})".format(var1gen=var1_config["gen"],var1reco=var1_config["reco"],var1gen_low=var1_config["binedgesgen"][i],var1gen_high=var1_config["binedgesgen"][i+1],var1reco_low=var1_config["binedgesreco"][j],var1reco_high=var1_config["binedgesreco"][j+1])
        mig[i][j].npy_cut=np_cuts[CutType.PassReco_PassGen]+[[var1_config["gen_key"],">=",str(var1_config["binedgesgen"][i])],[var1_config["gen_key"],"<",str(var1_config["binedgesgen"][i+1])],[var1_config["reco_key"],">=",str(var1_config["binedgesreco"][j])],[var1_config["reco_key"],"<",str(var1_config["binedgesreco"][j+1])]]
        if from_root:
          mig[i][j].fill_2Dhist_from_root(source,genWeight=genWeight)
        else:
          mig[i][j].fill_2Dhist_from_npz(files=source,weightarray=weight_array,genWeight=genWeight)
        mig[i][j].binwise_normalize_2Dhist()
        print "filled histogram:",mig[i][j].name
  else:
    mig = None

  hists["mig"] = mig
  return hists 

def get_sys_variations( config ):
    tree_sys_list = []
    fin_sys={}
    tree_sys={}
    for sys in config["inputfilesim_sys"].keys():
      fin_sys[sys]={
        "up":ROOT.TFile(config["inputfilesim_sys"][sys]["up"]),
        "down":ROOT.TFile(config["inputfilesim_sys"][sys]["down"]) if config["inputfilesim_sys"][sys]["down"] is not None else None
      }
      tree_sys[sys]={
        "up":fin_sys[sys]["up"].Get("ntuplizer/tree") if fin_sys[sys]["up"].Get("ntuplizer/tree") else fin_sys[sys]["up"].Get("tree"),
        "down":(fin_sys[sys]["down"].Get("ntuplizer/tree") if fin_sys[sys]["down"].Get("ntuplizer/tree") else fin_sys[sys]["down"].Get("tree") ) if fin_sys[sys]["down"] is not None else None
      }
    for sys in tree_sys.keys():
      tree_sys_list.append(tree_sys[sys]["up"])
      if tree_sys[sys]["down"] is not None:
        tree_sys_list.append(tree_sys[sys]["down"])

    return tree_sys_list

def get_bin_edges( conf, v1_dct, v2_dct, ttree, trees_syst):
    if conf['mergerecobin']:
      bin_edges_reco_merge=merge_bins(obs=[var1_dct["reco"],var2_dct["reco"]],trees=[ttree]+trees_syst,root_cut=root_cuts[CutType.PassReco_PassGen],threshold=conf["mergethresholdreco"],bin_edges_dim1_1d=var1_dct["binedgesreco"],bin_edges_dim2_1d=var2_dct["binedgesreco"])
    else:
      bin_edges_reco_merge=([var2_dct["binedgesreco"]]*var1_dct["nbinsreco"] if FineBin else [var2_dct["minreco"]+(var2_dct["maxreco"]-var2_dct["minreco"])/var2_dct["nbinsreco"]*ibinreco1  for ibinreco1 in range(var2_dct["nbinsreco"])]*var1_dct["nbinsreco"])

    if conf['mergegenbin']:
      bin_edges_gen_merge=merge_bins(obs=[var1_dct["gen"],var2_dct["gen"]],trees=[ttree]+trees_syst,root_cut=root_cuts[CutType.PassReco_PassGen],threshold=conf["mergethresholdgen"],bin_edges_dim1_1d=var1_dct["binedgesgen"],bin_edges_dim2_1d=var2_dct["binedgesgen"])
    else:
      bin_edges_gen_merge=([var2_dct["binedgesgen"]]*var1_dct["nbinsgen"] if FineBin else [var2_dct["mingen"]+(var2_dct["maxgen"]-var2_dct["mingen"])/var2_dct["nbinsgen"]*ibingen1  for ibingen1 in range(var2_dct["nbinsgen"])]*var1_dct["nbinsgen"])
    return bin_edges_reco_merge, bin_edges_gen_merge

def get_tree_data( config, pseudodata_NPZ):

    fin_refdata=ROOT.TFile(config["inputfiledata"],"READ")
    tree_refdata=fin_refdata.Get("ntuplizer/tree") if fin_refdata.Get("ntuplizer/tree") else fin_refdata.Get("tree")
    if config["pseudodata"]:
      if not pseudodata_NPZ:
        fin_data=ROOT.TFile(config["inputfilepseudodata"],"READ")
        tree_data = fin_data.Get("ntuplizer/tree") if  fin_data.Get("ntuplizer/tree") else  fin_data.Get("tree")
    else:
      tree_data = tree_refdata

    tree_data.SetDirectory(0)
    tree_refdata.SetDirectory(0)

    return tree_data, tree_refdata

if __name__=="__main__":

    parser = ArgumentParser()
    parser.add_argument('--config',default="Config_checkunfold/Config_sph_1d_v7_badtunesys.json",help="The configration file including the unfolding setup")
    parser.add_argument('--method',default="omnifold",help="omnifold/multfold/unifold")
    parser.add_argument('--migiter',type=int,nargs='+',help="Which iterations to plot the migration matrices")
    parser.add_argument('--step1',action="store_true",default=False,help="Process the histograms from the step1, otherwise process the step 2 results")
    parser.add_argument('--eff-acc',action="store_true",default=False,help="Consider the efficiency and acceptance in the unfolding")
    parser.add_argument('--eff-from-nominal',action="store_true",default=False,help="Apply the reconstruction efficiency of the nominal MC to the unfolded one, otherwise use the efficiency given by the unfolding algorithm.")
    args = parser.parse_args()

    with open(args.config, 'r') as configjson:
        config = json.load(configjson)
    with open(config["varunfold"], 'r') as fjson:
        info_var = json.load(fjson)

    var1_dct = info_var[config["var1"]]
    var2_dct = info_var[config["var2"]]

    if not os.path.exists(config["outputdir"]):
      os.makedirs(config["outputdir"])
    fin = ROOT.TFile(config["inputfilesim"],"READ")
    tree = fin.Get("ntuplizer/tree") if fin.Get("ntuplizer/tree") else fin.Get("tree")

    weightname=config["reweight"] if ("reweight" in config.keys() and config["reweight"]!="") else "genWeight"

    tree_sys_list=[]
    if config["addsys"]:
        tree_syst_list = get_sys_variations( config )

    FineBin = ("binedgesreco" in var1_dct.keys()) and ("binedgesreco" in var2_dct.keys())

    bin_edges_reco_merge, bin_edges_gen_merge = get_bin_edges( config, var1_dct, var2_dct, tree, tree_sys_list)

    mc_hists = fill_hist_lists("MC",var1_dct,var2_dct,bin_edges_gen_merge,bin_edges_reco_merge,tree,genWeight=weightname,from_root=True,weight_array=None,store_mig=True)

    #efficiency=reconstructed and generated / generated
    gen_inveff = hist_list("HistGenInvEff")
    gen_inveff.read_settings_from_config_dim1(var1_dct,isgen=True)
    gen_inveff.read_settings_from_config_dim2(var2_dct,isgen=True)
    gen_inveff.bin_edges_dim2 = bin_edges_gen_merge
    gen_inveff.fill_root_hists_name()
    gen_inveff.get_hist_from_division(mc_hists["gen_inclusive"],mc_hists["gen_passreco"])

    pseudodata_NPZ =  config["pseudodata"] and (isinstance(config["inputfilepseudodata"],list) or '.npz' in config["inputfilepseudodata"])
    tree_data, tree_refdata = get_tree_data( config, pseudodata_NPZ)

    if config["pseudodata"]:
      weight_pseudodata = None
      if pseudodata_NPZ:
        if config.get("pseudodataweight",None) is not None:
          file_weight_pseudodata=np.load(config["pseudodataweight"],allow_pickle=True)
          weight_pseudodata=file_weight_pseudodata[-1]

        event_data = config["inputfilepseudodata"]
        from_root = False
      else:
        event_data = tree_data
        from_root = True
      pseudo_hists = fill_hist_lists("Pseudodata",var1_dct,var2_dct,bin_edges_gen_merge,bin_edges_reco_merge,event_data,genWeight=weightname,from_root=from_root,weight_array=weight_pseudodata,store_mig=True)
      reco_data_tree = tree_refdata
      tag = "_Ref"
    else:
      reco_data_tree = tree_data
      tag = ""
    data_hists = fill_hist_lists("Data",var1_dct,var2_dct,None,bin_edges_reco_merge,reco_data_tree,genWeight="",from_root=True,weight_array=None,tag="_Ref", reco_only=True)

    if config["pseudodata"]:
      mc_norm_factor=pseudo_hists["reco_inclusive"].norm/mc_hists["reco_inclusive"].norm
    else:
      mc_norm_factor= data_hists["reco_inclusive"].norm/mc_hists["reco_inclusive"].norm

    for key, hist in mc_hists.items():
        if not (key == 'mig'):
            hist.multiply(mc_norm_factor)

    if not os.path.exists(config[args.method]["weight"]):
      print "Cannot find weight file ",config[args.method]["weight"]
      exit(0)

    weights=np.load(config[args.method]["weight"],allow_pickle=True)
    weights_per_iter = 4 if args.eff_acc else 2

    config["outputdir"]+"/unfold_"+config["var1"]+"_"+config["var2"]+"_"+config["MCtag"]+"_optimize_"+args.method+("_step1" if args.step1 else "")+".root"
    fout = ROOT.TFile(config["outputdir"]+"/unfold_"+config["var1"]+"_"+config["var2"]+"_"+config["MCtag"]+"_optimize_"+args.method+("_step1" if args.step1 else "")+".root","recreate")

    niter = len(weights)/weights_per_iter
    for i in range(0,niter+1):
      offset = weights_per_iter / 2 if (args.step1 and i > 0 ) else 0
      weight_iter = weights[weights_per_iter*i - offset ]

      store_mig = i in args.migiter

      unfold_hists = fill_hist_lists("MC_"+args.method,var1_dct,var2_dct,bin_edges_gen_merge,bin_edges_reco_merge,config[args.method]["sim"],genWeight=weightname,from_root=False,weight_array=weight_iter,store_mig=store_mig,tag="_iter"+str(i))
      if args.eff_from_nominal:
        gen_unfold_inclusive.get_hist_from_multiplication(unfold_hists["gen_passreco"],gen_inveff)
      if config["pseudodata"]:
        unf_norm_factor = pseudo_hists["reco_inclusive"].norm/unfold_hists["reco_inclusive"].norm
      else:
        unf_norm_factor = data_hists["reco_inclusive"].norm/unfold_hists["reco_inclusive"].norm

      for key, hist in unfold_hists.items():
        if key == 'mig':
            if i in args.migiter:
                for column in hist:
                    for row in column:
                        row.write_2Dhist()
        else:
            hist.multiply(unf_norm_factor)
            hist.write_hist_list()

    for key, hist in mc_hists.items():
        if key == 'mig':
          for column in hist:
            for row in column:
              row.write_2Dhist()

        else:
          hist.multiply(mc_norm_factor)
          hist.write_hist_list()

    gen_inveff.write_hist_list()
    if config["pseudodata"]:
      for key, hist in pseudo_hists.items():
        if key == 'mig':
            for column in hist:
              for row in column:
                row.write_2Dhist()
        else:
            hist.write_hist_list()
    else:
        data_hists["reco_inclusive"].write_hist_list()

