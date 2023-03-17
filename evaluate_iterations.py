import numpy as np
import json
from argparse import ArgumentParser
import os
import ROOT
from unfold_utils import *
from configs import ObsConfig

def fill_hist_lists(dataset,o1,o2,edges_gen,edges_reco,source,genWeight="",from_root=True,weight_array=None,store_mig=False,tag="", reco_only=False):
  hists = {}

  reco_inclusive=HistList("HistRecoInclusive_"+dataset,tag)
  print(o1.reco.edges)
  reco_inclusive.read_settings_from_config_dim1(o1,isgen=False)
  reco_inclusive.read_settings_from_config_dim2(o2,isgen=False)
  reco_inclusive.bin_edges_dim2 = edges_reco
  reco_inclusive.cut = cuts[CutType.PassReco]
  hists["reco_inclusive"] = reco_inclusive

  if not reco_only:
    gen_passreco=HistList("HistGen_"+dataset,tag)
    gen_passreco.read_settings_from_config_dim1(o1,isgen=True)
    gen_passreco.read_settings_from_config_dim2(o2,isgen=True)
    gen_passreco.bin_edges_dim2 = edges_gen
    gen_passreco.cut = cuts[CutType.PassReco_PassGen]
    hists["gen_passreco"] = gen_passreco

    gen_inclusive=HistList("HistGenInclusive_"+dataset,tag)
    gen_inclusive.read_settings_from_config_dim1(o1,isgen=True)
    gen_inclusive.read_settings_from_config_dim2(o2,isgen=True)
    gen_inclusive.bin_edges_dim2 = edges_gen
    gen_inclusive.cut = cuts[CutType.PassGen]
    hists["gen_inclusive"] = gen_inclusive

    reco_passgen=HistList("HistReco_"+dataset,tag)
    reco_passgen.read_settings_from_config_dim1(o1,isgen=False)
    reco_passgen.read_settings_from_config_dim2(o2,isgen=False)
    reco_passgen.bin_edges_dim2 = edges_reco
    reco_passgen.cut = cuts[CutType.PassReco_PassGen]
    hists["reco_passgen"] = reco_passgen

    gen_eff = HistList("HistGenEff_"+dataset,tag)
    gen_eff.read_settings_from_config_dim1(o1,isgen=True)
    gen_eff.read_settings_from_config_dim2(o2,isgen=True)
    gen_eff.bin_edges_dim2 = edges_gen
    gen_eff.fill_root_hists_name()

    reco_acc = HistList("HistRecoAcc_"+dataset,tag)
    reco_acc.read_settings_from_config_dim1(o1,isgen=False)
    reco_acc.read_settings_from_config_dim2(o2,isgen=False)
    reco_acc.bin_edges_dim2 = edges_reco
    reco_acc.fill_root_hists_name()

  for _, hist in hists.items():
    hist.fill_root_hists_name()
    print("histograms:", hist.root_hists_name)
    hist.fill_hist(source, from_root, weightarray=weight_array,genWeight=genWeight)
  print("filled")

  if not reco_only:
    gen_eff.get_hist_from_division(gen_passreco,gen_inclusive)
    reco_acc.get_hist_from_division(reco_passgen,reco_inclusive)
    hists["eff"] = gen_eff
    hists["acc"] = reco_acc

  mig = None
  if store_mig:
    mig = []
    for i, edge_i in enumerate(o1.gen.edges[:-1]):
      mig.append([])
      for j, edge_j in enumerate(o1.reco.edges[:-1]):
        edge_ip1 = o1.gen.edges[i+1]
        edge_jp1 = o1.reco.edges[j+1]

        #mig_ij
        this_mig = HistList(f'HistMig_{dataset}_{o1.gen.np_var}{edge_i}-{edge_ip1}_{o1.reco.np_var}{edge_j}-{edge_jp1}{tag}')
        mig[i].append(this_mig)

        mig[i][j].read_settings_from_config_dim1(o2,isgen=True)
        mig[i][j].read_settings_from_config_dim2(o2,isgen=False)

        gen_cuts_root = f" ({o1.gen.root_var}>={edge_i}) * ({o1.gen.root_var}<{edge_ip1}) "
        reco_cuts_root = f" ({o1.reco.root_var}>={edge_j}) * ({o1.reco.root_var}<{edge_jp1}) "
        all_cuts_root = root_cuts[CutType.PassReco_PassGen] + "*" + gen_cuts_root + "*" + reco_cuts_root
        mig[i][j].root_cut = all_cuts_root

        gen_cuts_np = [[o1.gen.np_var,">=",str(edge_i)], [o1.gen.np_var,"<",str(edge_ip1)]]
        reco_cuts_np = [[o1.reco.np_var, ">=", str(edge_j)], [o1.reco.np_var,"<",str(edge_jp1)]]
        all_cuts_np = np_cuts[CutType.PassReco_PassGen] + gen_cuts_np + reco_cuts_np
        mig[i][j].npy_cut = all_cuts_np

        mig[i][j].fill_2Dhist(source, from_root=from_root, weightarray=weight_array, genWeight=genWeight)
        mig[i][j].binwise_normalize_2Dhist()
        print("filled histogram:",mig[i][j].name)

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

def get_bin_edges( conf, var1, var2, trees):

    cut = root_cuts[CutType.PassReco_PassGen]
    bin_edges = []
    for lvl in ["reco","gen"]:
        do_merge = conf[f"merge{lvl}bin"]
        if do_merge:
            edges = merge_bins(obs=[var1[lvl].root_var , var2[lvl].root_var], trees=trees, root_cut=cut, threshold=conf[f"mergethreshold{lvl}"], bin_edges_dim1_1d=var1[lvl].edges, bin_edges_dim2_1d=var2[lvl].edges)
        else:
            if FineBin:
                edges = [var2[lvl].edges * var1[lvl].nbins ]
            else:
                edges = [ var2[lvl].min +(var2[lvl].max - var2[{lvl}].min)/ var2[lvl].nbins * ibin for ibin in range(var2[lvl].nbins)] * var1[lvl].nbins
        bin_edges.append(edges)

    return bin_edges



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

    obs1 = ObsConfig.from_yaml( config["varunfold"], [config["var1"]] )
    obs2 = ObsConfig.from_yaml( config["varunfold"], [config["var2"]] )
    #var1_dct = info_var[config["var1"]]
    #var2_dct = info_var[config["var2"]]
    print(obs1)
    print(obs2)

    if not os.path.exists(config["outputdir"]):
      os.makedirs(config["outputdir"])
    fin = ROOT.TFile(config["inputfilesim"],"READ")
    tree = fin.Get("ntuplizer/tree") if fin.Get("ntuplizer/tree") else fin.Get("tree")

    weightname=config["reweight"] if ("reweight" in config.keys() and config["reweight"]!="") else "genWeight"

    tree_sys_list=[]
    if config["addsys"]:
        tree_syst_list = get_sys_variations( config )

    #FineBin = ("binedgesreco" in var1_dct.keys()) and ("binedgesreco" in var2_dct.keys())
    FineBin = True

    all_trees = [tree] + tree_sys_list
    bin_edges_reco, bin_edges_gen = get_bin_edges( config, obs1, obs2, all_trees)

    mc_hists = fill_hist_lists("MC", obs1, obs2, bin_edges_gen, bin_edges_reco, tree, genWeight=weightname, store_mig=True)

    #efficiency=reconstructed and generated / generated
    gen_inveff = HistList("HistGenInvEff")
    gen_inveff.read_settings_from_config_dim1(obs1,isgen=True)
    gen_inveff.read_settings_from_config_dim2(obs2,isgen=True)
    gen_inveff.bin_edges_dim2 = bin_edges_gen
    gen_inveff.fill_root_hists_name()
    gen_inveff.get_hist_from_division(mc_hists["gen_inclusive"],mc_hists["gen_passreco"])

    pseudodata_NPZ =  config["pseudodata"] and (isinstance(config["inputfilepseudodata"],list) or '.npz' in config["inputfilepseudodata"])
    #tree_data, tree_refdata = get_tree_data( config, pseudodata_NPZ)

    tree_data = None
    fin_refdata=ROOT.TFile(config["inputfiledata"],"READ")
    tree_refdata=fin_refdata.Get("ntuplizer/tree") if fin_refdata.Get("ntuplizer/tree") else fin_refdata.Get("tree")
    if config["pseudodata"]:
      if not pseudodata_NPZ:
        fin_data=ROOT.TFile(config["inputfilepseudodata"],"READ")
        tree_data = fin_data.Get("ntuplizer/tree") if  fin_data.Get("ntuplizer/tree") else  fin_data.Get("tree")
    else:
      tree_data = tree_refdata

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
      pseudo_hists = fill_hist_lists("Pseudodata", obs1, obs2, bin_edges_gen, bin_edges_reco, event_data, genWeight=weightname, from_root=from_root, weight_array=weight_pseudodata,store_mig=True)
      reco_data_tree = tree_refdata
      tag = "_Ref"
    else:
      reco_data_tree = tree_data
      tag = ""
    data_hists = fill_hist_lists("Data", obs1, obs2, None, bin_edges_reco, reco_data_tree, tag="_Ref", reco_only=True)

    normalization_hist = pseudo_hists["reco_inclusive"] if config["pseudodata"] else data_hists["reco_inclusive"]
    mc_norm_factor = normalization_hist.norm / mc_hists["reco_inclusive"].norm

    for key, hist in mc_hists.items():
        if not (key == 'mig' or key == "eff" or key == "acc"):
            hist.multiply(mc_norm_factor)

    if not os.path.exists(config[args.method]["weight"]):
      print("Cannot find weight file ",config[args.method]["weight"])
      exit(0)

    weights=np.load(config[args.method]["weight"],allow_pickle=True)
    weights_per_iter = 4 if args.eff_acc else 2

    step1_tag = "_step1" if args.step1 else ""
    fout = ROOT.TFile(f'{config["outputdir"]}/unfold_{config["var1"]}_{config["var2"]}_{config["MCtag"]}_optimize_{args.method}{step1_tag}.root', "recreate")

    niter = len(weights)//weights_per_iter
    for i in range(0,niter+1):
      offset = weights_per_iter // 2 if (args.step1 and i > 0 ) else 0
      weight_iter = weights[weights_per_iter*i - offset ]
      store_mig = i in args.migiter

      unfold_hists = fill_hist_lists("MC_"+args.method, obs1, obs2, bin_edges_gen,bin_edges_reco,config[args.method]["sim"],genWeight=weightname,from_root=False,weight_array=weight_iter,store_mig=store_mig,tag="_iter"+str(i))
      if args.eff_from_nominal:
        unfold_hists["gen_inclusive"].get_hist_from_multiplication(unfold_hists["gen_passreco"],gen_inveff)
      unf_norm_factor = normalization_hist.norm / unfold_hists["reco_inclusive"].norm

      for key, hist in unfold_hists.items():
        if not (key == 'mig' or key == "eff" or key == "acc"):
            hist.multiply(unf_norm_factor)

      write_all_hists(unfold_hists)

    write_all_hists(mc_hists)
    gen_inveff.write_hist_list()
    if config["pseudodata"]:
      write_all_hists(pseudo_hists)
    else:
      data_hists["reco_inclusive"].write_hist_list()

