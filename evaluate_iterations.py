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

  mig = None
  if store_mig:
    mig = []
    for i, edge_i in enumerate(var1_config["binedgesgen"][:-1]):
      mig.append([])
      for j, edge_j in enumerate(var1_config["binedgesreco"][:-1]):
        edge_ip1 = var1_config["binedgesgen"][i+1]
        edge_jp1 = var1_config["binedgesreco"][j+1]

        #mig_ij
        this_mig = hist_list("HistMig_" +dataset+"_" + var1_config["gen_key"] +str(edge_i) + "-" +str(edge_ip1) + "_" +var1_config["reco_key"] + str(edge_j) + "-" + str(edge_jp1) + tag )
        mig[i].append(this_mig)

        mig[i][j].read_settings_from_config_dim1(var2_config,isgen=True)
        mig[i][j].read_settings_from_config_dim2(var2_config,isgen=False)

        gen_cuts_root = " ({}>={}) * ({}<{}) ".format(var1_config["gen"], edge_i, var1_config["gen"], edge_ip1)
        reco_cuts_root = " ({}>={}) * ({}<{}) ".format(var1_config["reco"], edge_j, var1_config["reco"], edge_jp1)
        all_cuts_root = root_cuts[CutType.PassReco_PassGen] + "*" + gen_cuts_root + "*" + reco_cuts_root
        mig[i][j].root_cut = all_cuts_root

        gen_cuts_np = [[var1_config["gen_key"],">=",str(edge_i)], [var1_config["gen_key"],"<",str(edge_ip1)]]
        reco_cuts_np = [[var1_config["reco_key"], ">=", str(edge_j)], [var1_config["reco_key"],"<",str(edge_jp1)]]
        all_cuts_np = np_cuts[CutType.PassReco_PassGen] + gen_cuts_np + reco_cuts_np
        mig[i][j].npy_cut = all_cuts_np

        mig[i][j].fill_2Dhist(source, from_root=from_root, weightarray=weight_array, genWeight=genWeight)
        mig[i][j].binwise_normalize_2Dhist()
        print "filled histogram:",mig[i][j].name

  hists["mig"] = mig
  return hists

def get_sys_variations( config ):
    tree_sys_list = []
    for sys in config["inputfilesim_sys"].keys():
      for vari in ["up", "down"]:
        fin_sys_vari = ROOT.TFile(config["inputfilesim_sys"][sys][vari])
        default = fin_sys_vari.Get("ntuplizer/tree")
        if default:
            tree = default
        else:
            tree = fin_sys_vari.Get("tree")
        if tree is not None:
            tree_sys_list.append(tree)
    return tree_sys_list

def get_bin_edges( conf, v1_dct, v2_dct, trees):

    cut = root_cuts[CutType.PassReco_PassGen]
    bin_edges = []
    for lvl in ["reco","gen"]:
        do_merge = conf["merge" + lvl + "bin"]
        if do_merge:
            edges = merge_bins(obs=[v1_dct[lvl] , v2_dct[lvl]], trees=trees, root_cut=cut, threshold=conf["mergethreshold" + lvl], bin_edges_dim1_1d=v1_dct["binedges" +lvl], bin_edges_dim2_1d=v2_dct["binedges" + lvl])
        else:
            if FineBin:
                edges = [v2_dct["binedges" + lvl] * v1_dct["nbins" + lvl] ]
            else:
                edges = [ v2_dct["min"+lvl]+(v2_dct["max" + lvl] - v2_dct["min"+lvl])/ v2_dct["nbins" + lvl] * ibin for ibin in range(v2_dct["nbins" + lvl])] * v1_dct["nbins" + lvl]
        bin_edges.append(edges)

    return bin_edges

def get_tree_data( config, pseudodata_NPZ):

    fin_refdata=ROOT.TFile(config["inputfiledata"],"READ")
    tree_data = None
    tree_refdata=fin_refdata.Get("ntuplizer/tree") if fin_refdata.Get("ntuplizer/tree") else fin_refdata.Get("tree")
    if config["pseudodata"]:
      if not pseudodata_NPZ:
        fin_data=ROOT.TFile(config["inputfilepseudodata"],"READ")
        tree_data = fin_data.Get("ntuplizer/tree") if  fin_data.Get("ntuplizer/tree") else  fin_data.Get("tree")
        tree_data.SetDirectory(0)
    else:
      tree_data = tree_refdata
      tree_data.SetDirectory(0)

    tree_refdata.SetDirectory(0)

    return tree_data, tree_refdata

def write_all_hists( hist_dict ):
    for key, hist in hist_dict.items():
        if hist == None:
            continue

        if key == 'mig':
          for column in hist:
            for row in column:
              row.write_2Dhist()
        else:
          hist.write_hist_list()

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

    all_trees = [tree] + tree_sys_list
    bin_edges_reco, bin_edges_gen = get_bin_edges( config, var1_dct, var2_dct, all_trees)

    mc_hists = fill_hist_lists("MC", var1_dct, var2_dct, bin_edges_gen, bin_edges_reco, tree, genWeight=weightname, store_mig=True)

    #efficiency=reconstructed and generated / generated
    gen_inveff = hist_list("HistGenInvEff")
    gen_inveff.read_settings_from_config_dim1(var1_dct,isgen=True)
    gen_inveff.read_settings_from_config_dim2(var2_dct,isgen=True)
    gen_inveff.bin_edges_dim2 = bin_edges_gen
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
      pseudo_hists = fill_hist_lists("Pseudodata",var1_dct,var2_dct,bin_edges_gen,bin_edges_reco,event_data,genWeight=weightname,from_root=from_root,weight_array=weight_pseudodata,store_mig=True)
      reco_data_tree = tree_refdata
      tag = "_Ref"
    else:
      reco_data_tree = tree_data
      tag = ""
    data_hists = fill_hist_lists("Data", var1_dct, var2_dct, None, bin_edges_reco, reco_data_tree, tag="_Ref", reco_only=True)

    normalization_hist = pseudo_hists["reco_inclusive"] if config["pseudodata"] else data_hists["reco_inclusive"]
    mc_norm_factor = normalization_hist.norm / mc_hists["reco_inclusive"].norm

    for key, hist in mc_hists.items():
        if not (key == 'mig'):
            hist.multiply(mc_norm_factor)

    if not os.path.exists(config[args.method]["weight"]):
      print "Cannot find weight file ",config[args.method]["weight"]
      exit(0)

    weights=np.load(config[args.method]["weight"],allow_pickle=True)
    weights_per_iter = 4 if args.eff_acc else 2

    fout = ROOT.TFile(config["outputdir"]+"/unfold_"+config["var1"]+"_"+config["var2"]+"_"+config["MCtag"]+"_optimize_"+args.method+("_step1" if args.step1 else "")+".root","recreate")

    niter = len(weights)/weights_per_iter
    for i in range(0,niter+1):
      offset = weights_per_iter / 2 if (args.step1 and i > 0 ) else 0
      weight_iter = weights[weights_per_iter*i - offset ]
      store_mig = i in args.migiter

      unfold_hists = fill_hist_lists("MC_"+args.method,var1_dct,var2_dct,bin_edges_gen,bin_edges_reco,config[args.method]["sim"],genWeight=weightname,from_root=False,weight_array=weight_iter,store_mig=store_mig,tag="_iter"+str(i))
      if args.eff_from_nominal:
        gen_unfold_inclusive.get_hist_from_multiplication(unfold_hists["gen_passreco"],gen_inveff)
      unf_norm_factor = normalization_hist.norm / unfold_hists["reco_inclusive"].norm

      write_all_hists(unfold_hists)

    write_all_hists(mc_hists)
    gen_inveff.write_hist_list()
    if config["pseudodata"]:
      write_all_hists(pseudo_hists)
    else:
      data_hists["reco_inclusive"].write_hist_list()

