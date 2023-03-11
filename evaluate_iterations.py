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

def fill_hist_lists(dataset,var1_config,var2_config,edges_gen,edges_reco,source,genWeight="",from_root=True,weight_array=None,store_mig=False,tag=""):
  gen_passreco=hist_list("HistGen_"+dataset,tag)
  gen_passreco.read_settings_from_config_dim1(var1_config,isgen=True)
  gen_passreco.read_settings_from_config_dim2(var2_config,isgen=True)
  gen_passreco.bin_edges_dim2 = edges_gen
  gen_passreco.cut = cuts[CutType.PassReco_PassGen]
  gen_passreco.fill_root_hists_name()
  print "histograms:",gen_passreco.root_hists_name

  gen_inclusive=hist_list("HistGenInclusive_"+dataset,tag)
  gen_inclusive.read_settings_from_config_dim1(var1_config,isgen=True)
  gen_inclusive.read_settings_from_config_dim2(var2_config,isgen=True)
  gen_inclusive.bin_edges_dim2 = edges_gen
  gen_inclusive.cut = cuts[CutType.PassGen]
  gen_inclusive.fill_root_hists_name()
  print "histograms:",gen_inclusive.root_hists_name

  reco_passgen=hist_list("HistReco_"+dataset,tag)
  reco_passgen.read_settings_from_config_dim1(var1_config,isgen=False)
  reco_passgen.read_settings_from_config_dim2(var2_config,isgen=False)
  reco_passgen.bin_edges_dim2 = edges_reco
  reco_passgen.cut = cuts[CutType.PassReco_PassGen]
  reco_passgen.fill_root_hists_name()
  print "histograms:",reco_passgen.root_hists_name

  reco_inclusive=hist_list("HistRecoInclusive_"+dataset,tag)
  reco_inclusive.read_settings_from_config_dim1(var1_config,isgen=False)
  reco_inclusive.read_settings_from_config_dim2(var2_config,isgen=False)
  reco_inclusive.bin_edges_dim2 = edges_reco
  reco_inclusive.cut = cuts[CutType.PassReco]
  reco_inclusive.fill_root_hists_name()
  print "histograms:",reco_inclusive.root_hists_name

  if from_root:
    gen_passreco.fill_hist_from_root(source,genWeight=genWeight)
    gen_inclusive.fill_hist_from_root(source,genWeight=genWeight)
    reco_passgen.fill_hist_from_root(source,genWeight=genWeight)
    reco_inclusive.fill_hist_from_root(source,genWeight=genWeight)
  else:
    gen_passreco.fill_hist_from_npz(files=source,weightarray=weight_array,genWeight=genWeight)
    gen_inclusive.fill_hist_from_npz(files=source,weightarray=weight_array,genWeight=genWeight)
    reco_passgen.fill_hist_from_npz(files=source,weightarray=weight_array,genWeight=genWeight)
    reco_inclusive.fill_hist_from_npz(files=source,weightarray=weight_array,genWeight=genWeight)
  print "filled"
  if store_mig:
    mig = [[hist_list("HistMig_"+dataset+"_"+var1_config["gen_key"]+str(var1_config["binedgesgen"][i])+"-"+str(var1_config["binedgesgen"][i+1])+"_"+var1_config["reco_key"]+str(var1_config["binedgesreco"][j])+"-"+str(var1_config["binedgesreco"][j+1])+tag) for j in range(len(var1_config["binedgesreco"])-1)] for i in range(len(var1_config["binedgesgen"])-1)]
    for i in range(len(var1_config["binedgesgen"])-1):
      for j in range(len(var1_config["binedgesreco"])-1):
        mig[i][j].read_settings_from_config_dim1(var2_config,isgen=True)
        mig[i][j].read_settings_from_config_dim2(var2_config,isgen=False)
        mig[i][j].root_cut=root_cut_passreco_passgen+"*({var1gen}>={var1gen_low})*({var1gen}<{var1gen_high})*({var1reco}>={var1reco_low})*({var1reco}<{var1reco_high})".format(var1gen=var1_config["gen"],var1reco=var1_config["reco"],var1gen_low=var1_config["binedgesgen"][i],var1gen_high=var1_config["binedgesgen"][i+1],var1reco_low=var1_config["binedgesreco"][j],var1reco_high=var1_config["binedgesreco"][j+1])
        mig[i][j].npy_cut=np_cut_passreco_passgen+[[var1_config["gen_key"],">=",str(var1_config["binedgesgen"][i])],[var1_config["gen_key"],"<",str(var1_config["binedgesgen"][i+1])],[var1_config["reco_key"],">=",str(var1_config["binedgesreco"][j])],[var1_config["reco_key"],"<",str(var1_config["binedgesreco"][j+1])]]
        if from_root:
          mig[i][j].fill_2Dhist_from_root(source,genWeight=genWeight)
        else:
          mig[i][j].fill_2Dhist_from_npz(files=source,weightarray=weight_array,genWeight=genWeight)
        mig[i][j].binwise_normalize_2Dhist()
        print "filled histogram:",mig[i][j].name
  else:
    mig = None
  return gen_passreco,gen_inclusive,reco_passgen,reco_inclusive,mig

def fill_hist_lists_recoonly(dataset,var1_config,var2_config,edges_reco,source,genWeight="",from_root=True,weight_array=None,tag=""):
  reco_inclusive=hist_list("HistRecoInclusive_"+dataset,tag)
  reco_inclusive.read_settings_from_config_dim1(var1_config,isgen=False)
  reco_inclusive.read_settings_from_config_dim2(var2_config,isgen=False)
  reco_inclusive.bin_edges_dim2 = edges_reco
  reco_inclusive.cut = cuts[CutType.PassReco]
  reco_inclusive.fill_root_hists_name()
  if from_root:
    reco_inclusive.fill_hist_from_root(source,genWeight=genWeight)
  else:
    reco_inclusive.fill_hist_from_npz(files=source,weightarray=weight_array,genWeight=genWeight)
  print "filled histogram:",reco_inclusive.root_hists_name
  return reco_inclusive


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
    if config['mergerecobin']:
      bin_edges_reco_merge=merge_bins(obs=[var1_dct["reco"],var2_dct["reco"]],trees=[ttree]+trees_syst,root_cut=root_cut_passreco_passgen,threshold=config["mergethresholdreco"],bin_edges_dim1_1d=var1_dct["binedgesreco"],bin_edges_dim2_1d=var2_dct["binedgesreco"])
    else:
      bin_edges_reco_merge=([var2_dct["binedgesreco"]]*var1_dct["nbinsreco"] if FineBin else [var2_dct["minreco"]+(var2_dct["maxreco"]-var2_dct["minreco"])/var2_dct["nbinsreco"]*ibinreco1  for ibinreco1 in range(var2_dct["nbinsreco"])]*var1_dct["nbinsreco"])

    if config['mergegenbin']:
      bin_edges_gen_merge=merge_bins(obs=[var1_dct["gen"],var2_dct["gen"]],trees=[ttree]+trees_syst,root_cut=root_cut_passreco_passgen,threshold=config["mergethresholdgen"],bin_edges_dim1_1d=var1_dct["binedgesgen"],bin_edges_dim2_1d=var2_dct["binedgesgen"])
    else:
      bin_edges_gen_merge=([var2_dct["binedgesgen"]]*var1_dct["nbinsgen"] if FineBin else [var2_dct["mingen"]+(var2_dct["maxgen"]-var2_dct["mingen"])/var2_dct["nbinsgen"]*ibingen1  for ibingen1 in range(var2_dct["nbinsgen"])]*var1_dct["nbinsgen"])
    return bin_edges_reco_merge, bin_edges_gen_merge

def get_tree_data( config, pseudodata_NPZ):
    if config["pseudodata"]:
      if not pseudodata_NPZ:
        fin_data=ROOT.TFile(config["inputfilepseudodata"],"READ")
        tree_data = fin_data.Get("ntuplizer/tree") if  fin_data.Get("ntuplizer/tree") else  fin_data.Get("tree")
      fin_refdata=ROOT.TFile(config["inputfiledata"],"READ")
      tree_refdata=fin_refdata.Get("ntuplizer/tree") if fin_refdata.Get("ntuplizer/tree") else fin_refdata.Get("tree")
    else:
      fin_data = ROOT.TFile(config["inputfiledata"],"READ")
      tree_data = fin_data.Get("ntuplizer/tree") if fin_data.Get("ntuplizer/tree") else fin_data.Get("tree")
      tree_refdata = None

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

    gen_MC_passreco,gen_MC_inclusive,reco_MC_passgen,reco_MC_inclusive,mig_MC=fill_hist_lists("MC",var1_dct,var2_dct,bin_edges_gen_merge,bin_edges_reco_merge,tree,genWeight=weightname,from_root=True,weight_array=None,store_mig=True)

    #efficiency=reconstructed and generated / generated
    gen_inveff = hist_list("HistGenInvEff")
    gen_inveff.read_settings_from_config_dim1(var1_dct,isgen=True)
    gen_inveff.read_settings_from_config_dim2(var2_dct,isgen=True)
    gen_inveff.bin_edges_dim2 = bin_edges_gen_merge
    gen_inveff.fill_root_hists_name()
    gen_inveff.get_hist_from_division(gen_MC_inclusive,gen_MC_passreco)

    pseudodata_NPZ =  config["pseudodata"] and (isinstance(config["inputfilepseudodata"],list) or '.npz' in config["inputfilepseudodata"])
    tree_data, tree_refdata = get_tree_data( config, pseudodata_NPZ)

    if config["pseudodata"]:
      if pseudodata_NPZ:
        if "pseudodataweight" in config.keys() and config["pseudodataweight"] is not None:
          file_weight_pseudodata=np.load(config["pseudodataweight"],allow_pickle=True)
          weight_pseudodata=file_weight_pseudodata[-1]
        else:
          weight_pseudodata=None
        gen_pseudodata_passreco,gen_pseudodata_inclusive,reco_pseudodata_passgen,reco_pseudodata_inclusive,mig_pseudodata=fill_hist_lists("Pseudodata",var1_dct,var2_dct,bin_edges_gen_merge,bin_edges_reco_merge,config["inputfilepseudodata"],genWeight=weightname,from_root=False,weight_array=weight_pseudodata,store_mig=True)
      else:
        gen_pseudodata_passreco,gen_pseudodata_inclusive,reco_pseudodata_passgen,reco_pseudodata_inclusive,mig_pseudodata=fill_hist_lists("Pseudodata",var1_dct,var2_dct,bin_edges_gen_merge,bin_edges_reco_merge,tree_data,genWeight=weightname,from_root=True,weight_array=None,store_mig=True)
      reco_data_inclusive=fill_hist_lists_recoonly("Data",var1_dct,var2_dct,bin_edges_reco_merge,tree_refdata,genWeight="",from_root=True,weight_array=None,tag="_Ref")
    else:
      reco_data_inclusive=fill_hist_lists_recoonly("Data",var1_dct,var2_dct,bin_edges_reco_merge,tree_data,genWeight="",from_root=True,weight_array=None)

    if config["pseudodata"]:
      norm_factor=reco_pseudodata_inclusive.norm/reco_MC_inclusive.norm
    else:
      norm_factor=reco_data_inclusive.norm/reco_MC_inclusive.norm

    gen_MC_passreco.multiply(norm_factor)
    gen_MC_inclusive.multiply(norm_factor)
    reco_MC_passgen.multiply(norm_factor)
    reco_MC_inclusive.multiply(norm_factor)

    if not os.path.exists(config[args.method]["weight"]):
      print "Cannot find weight file ",config[args.method]["weight"]
      exit(0)

    weights=np.load(config[args.method]["weight"],allow_pickle=True)
    if args.eff_acc:
      niter=len(weights)/4
    else:
      niter=len(weights)/2

    fout = ROOT.TFile(config["outputdir"]+"/unfold_"+config["var1"]+"_"+config["var2"]+"_"+config["MCtag"]+"_optimize_"+args.method+("_step1" if args.step1 else "")+".root","recreate")

    for i in range(0,niter+1):
      if args.eff_acc:
        if args.step1 and i>0:
          weight_iter=weights[4*i-2]
        else:
          weight_iter=weights[4*i]
      else:
        if args.step1 and i>0:
          weight_iter=weights[2*i-1]
        else:
          weight_iter=weights[2*i]
      gen_unfold_passreco,gen_unfold_inclusive,reco_unfold_passgen,reco_unfold_inclusive,mig_unfold=fill_hist_lists("MC_"+args.method,var1_dct,var2_dct,bin_edges_gen_merge,bin_edges_reco_merge,config[args.method]["sim"],genWeight=weightname,from_root=False,weight_array=weight_iter,store_mig=True if i in args.migiter else False,tag="_iter"+str(i))
      if args.eff_from_nominal:
        gen_unfold_inclusive.get_hist_from_multiplication(gen_unfold_passreco,gen_inveff)
      if config["pseudodata"]:
        norm_factor=reco_pseudodata_inclusive.norm/reco_unfold_inclusive.norm
      else:
        norm_factor=reco_data_inclusive.norm/reco_unfold_inclusive.norm
      gen_unfold_passreco.multiply(norm_factor)
      gen_unfold_inclusive.multiply(norm_factor)
      reco_unfold_passgen.multiply(norm_factor)
      reco_unfold_inclusive.multiply(norm_factor)
      gen_unfold_passreco.write_hist_list()
      gen_unfold_inclusive.write_hist_list()
      reco_unfold_passgen.write_hist_list()
      reco_unfold_inclusive.write_hist_list()
      if i in args.migiter:
        for column in mig_unfold:
          for row in column:
            row.write_2Dhist()
    gen_MC_passreco.write_hist_list()
    gen_MC_inclusive.write_hist_list()
    reco_MC_passgen.write_hist_list()
    reco_MC_inclusive.write_hist_list()
    for column in mig_MC:
      for row in column:
        row.write_2Dhist()
    gen_inveff.write_hist_list()
    if config["pseudodata"]:
      gen_pseudodata_passreco.write_hist_list()
      gen_pseudodata_inclusive.write_hist_list()
      reco_pseudodata_passgen.write_hist_list()
      reco_pseudodata_inclusive.write_hist_list()
      for column in mig_pseudodata:
        for row in column:
          row.write_2Dhist()
    else:
      reco_data_inclusive.write_hist_list()





