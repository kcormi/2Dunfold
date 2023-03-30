import numpy as np
import json
import yaml
from argparse import ArgumentParser
import os
import ROOT as rt
from Plotting_cfg import *
from scipy.stats import chi2 as CHI2
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from plot_utils import *
from unfold_utils import *
from arg_parsing import *
from plot_mplhep_utils import *
from dataclasses import dataclass, field, asdict
from typing import Union

from plot_configs import ResultPlotSettings, HistConfig, PlotConfig
from configs import ObsConfig
import pandas as pd

@dataclass
class HistConfigCollection:
  gen: HistConfig = None
  reco: HistConfig = None
  eff: HistConfig = None
  acc: HistConfig = None
  mig: HistConfig = None


def get_methods(config):
  methods = []
  if ("addomnifold" in config.keys()) and config["addomnifold"]:
    methods.append("omnifold")
  if ("addmultifold" in config.keys()) and config["addmultifold"]:
    methods.append("multifold")
  if ("addunifold" in config.keys()) and config["addunifold"]:
    methods.append("unifold")
  if ("addcombine" in config.keys()) and config["addcombine"]:
    methods.append("MLE")
  return methods

def get_df_entries(df, **kwargs_selections):
  for ikey,key in enumerate(kwargs_selections.keys()):
    if ikey == 0:
      selection = (df[key]==kwargs_selections[key])
    else:
      selection = selection & (df[key]==kwargs_selections[key])
    
  return df[selection]

def get_histconfigcollection(df,color_list,style_list,legend_list):
  gen_entry = df[(df['histtype']=='inclusive')&df['dim1_isgen']&df['dim2_isgen']]
  reco_entry = df[(df['histtype']=='inclusive')&~df['dim1_isgen']&~df['dim2_isgen']] 
  #Plan: process the results from bootstraps -> add 'histtype' options 'inclusive_sysunc', 'inclusive_statunc','inclusive_totalunc'
  eff_entry = df[(df['histtype']=='efficiency')&df['dim1_isgen']&df['dim2_isgen']]
  acc_entry = df[(df['histtype']=='acceptance')&~df['dim1_isgen']&~df['dim2_isgen']]
  mig_entry = df[(df['histtype']=='migration')&df['dim1_isgen']&~df['dim2_isgen']] 

  #print(gen_entry,reco_entry,eff_entry,acc_entry,mig_entry)
  if not isinstance(color_list,list):
    color_list = [color_list]*4+[None]

  if not isinstance(style_list,list):
    style_list = [style_list]*4+[None]

  if not isinstance(legend_list,list):
    legend_list = [legend_list,legend_list,legend_list+" Eff.",legend_list+" Eff.",legend_list+" Mig."]

  histarray_statonly = [0]*5 # Replace with the HistArray from stat-unc-only unfolding results

  list_HistConfig = []
  for i,entry in enumerate([gen_entry,reco_entry,eff_entry,acc_entry,mig_entry]):
    if len(entry)>0:
      histarray = HistArray.from_df_entry(entry.iloc[0])
      if i<2:
        histarray.divide_by_bin_width()
      list_HistConfig.append(HistConfig(histarray,
                                        histarray_statonly[i],
                                        color_list[i], style_list[i],legend_list[i]))
    else:
      list_HistConfig.append(None)

  return HistConfigCollection(*list_HistConfig)

def plot_wrapper(plt_list,**kwargs):
  input_args = {
    "hist_ref_stat": plt_list.ref.stat,
    "style_ref": plt_list.ref.style,
    "color_ref": plt_list.ref.color,
    "list_style_compare": plt_list.styles,
    "list_color_compare": plt_list.colors,
    "label_ratio":plt_list.ratio
  }
  use_root  = kwargs.pop("use_root")
  input_args.update(kwargs)
  if use_root:
    plot_flat_hists_root(plt_list.ref.hist, plt_list.hists, plt_list.ref.legend, plt_list.legends, **input_args)
  else:
    plot_flat_hists_mpl(plt_list.ref.hist, plt_list.hists, plt_list.ref.legend, plt_list.legends, **input_args)

def draw_plot(plt_list, plotdir, var1_nm, var2_nm, obs2, txt_list,use_root=True):
    if plt_list.ref is None: return
    path = f'{plotdir}/{plt_list.name}_{var1_nm}_{var2_nm}'
    os.system(f"mkdir -p {plotdir}")
    axis_title = obs2.reco.name if plt_list.is_reco else obs2.gen.name

    plot_args = {
      "title": axis_title,
      "is_logY": 0,
      "do_ratio": 1,
      "output_path": path,
      "text_list": txt_list,
      "labelY":'Normalized Events/Bin Width',
      "range_ratio": 0.1 if "eff" in plt_list.name or "acc" in plt_list.name else 0.3,
      "use_root": use_root
    }
    plot_wrapper(plt_list,**plot_args)

    if not("eff" in plt_list.name or "acc" in plt_list.name):
      plot_args["is_logY"] = 1
      plot_wrapper(plt_list,**plot_args)

if __name__=="__main__":

    parser = ArgumentParser()
    parser.add_argument('--config',default="config/plot_1d_v7_MCEPOS_unfoldCP1.json",help="The configration file including the unfolding setup")
    parser.add_argument('--config-style', type=str, default = 'config/results_style.yml')
    parser.add_argument('--obs', type=obs_pair,nargs = '+', help="Pair of observables separated by a comma")
    args = parser.parse_args()

    with open(args.config, 'r') as configjson:
        config = json.load(configjson)

    style_file = args.config_style
    with open(style_file, 'r') as style_json:
        config_style = yaml.safe_load(style_json)

    df = pd.read_csv(config["input"])



    for method in get_methods(config):
      chi2_records_dict = {}
      result_settings = ResultPlotSettings.from_yaml(style_file, ['mpl', method])
      result_settings.sys_reweight = False if config['workflow']=="unfold" else True
      for obs in args.obs:
      
        (obs1_name, obs1_binning), (obs2_name, obs2_binning) = parse_obs_args( config, obs )

        obs1 = ObsConfig.from_yaml( config["varunfold"], [obs1_name], binning=obs1_binning )
        obs2 = ObsConfig.from_yaml( config["varunfold"], [obs2_name], binning=obs2_binning )
        TextListReco = [f'{obs1.reco.edges[i]} #leq {obs1.reco.shortname} < {obs1.reco.edges[i+1]}' for i in range(obs1.reco.nbins)]
        TextListGen = [f'{obs1.gen.edges[i]} #leq {obs1.gen.shortname} < {obs1.gen.edges[i+1]}' for i in range(obs1.gen.nbins)]+(["background"] if config["addbkg"] else [])

        chi2_records_dict[f'{obs[0]}_{obs[1]}'] = {}
        base_sel = {'obs1': obs[0],
                    'obs2': obs[1],
                    'method': method
                   }
        gen_sel = {'dim1_isgen': True, 'dim2_isgen': True}
        reco_sel = {'dim1_isgen': False, 'dim2_isgen': False}
        mig_sel = {'dim1_isgen': True, 'dim2_isgen': False, 'histtype': 'migration' }
        if "pseudodata" in config.keys() and config["pseudodata"]:
          data_sel = {'datatype': "pseudodata",
                      'dataset': config["pseudodata_name"]}
          data_legend = config["pseudodatalegend"]
          data_color = config_style['mpl']["pseudodata_color"]
        else:
          data_sel = {'datatype': "data",
                      'dataset': config["data_name"]}
          data_legend = "Data"
          data_color = config_style['mpl']["data_color"]

        MC_sel = {'datatype': "MC",
                  'dataset': config["MC_name"]}
        MC_legend = config["MClegend"]
        MC_color = config_style['mpl']["MC_color"]

        records_data = get_df_entries(df, **base_sel|data_sel)
        records_MC = get_df_entries(df, **base_sel|MC_sel)

        hcc_data = get_histconfigcollection(records_data,data_color,['triangle','cross','triangle','triangle',None],data_legend)

        hcc_MC = get_histconfigcollection(records_MC,MC_color,['fillederror','fillederror','cross','cross',None],MC_legend)

        iters = np.sort(np.unique(df[~np.isnan(df['iter'])]['iter']))

        for it in iters:
          iter_sel = {'datatype': config['workflow'],
                      'dataset': config['MC_name'],
                      'iter': it}
          records_it = get_df_entries(df, **base_sel|iter_sel)
          it_color = result_settings.color
          it_color_unfold = result_settings.color_unfold
          it_legend = result_settings.legend
          it_legend_unfold = result_settings.legend_unfold
          tag = result_settings.tag
          hcc_it = get_histconfigcollection(records_it,[it_color_unfold,it_color,it_color_unfold,it_color_unfold,None],["marker","filled","cross","cross",None],[it_legend_unfold,it_legend,it_legend_unfold+" Eff.",it_legend_unfold+" Acc.",it_legend_unfold+" Mig."])

          pltLists = {}
          ratio = 'Refold / data' if config['workflow']=="unfold" else "Reweight / sys. var."
          pltLists[f"Refoldcompare_iter{it}"] = PlotConfig(
                                                           hcc_data.reco,
                                                           [hcc_it.reco,hcc_MC.reco],
                                                           f'data_refold_{tag}_iter{it}',
                                                           ratio)

          ratio = 'Unfold / MC' if config['workflow']=="unfold" else "Reweight / MC"
          pltLists[f"Unfoldcompare_iter{it}"] = PlotConfig(
                                                           hcc_MC.gen,
                                                           [hcc_it.gen],
                                                           f'MC_unfold_{tag}_iter{it}',
                                                           ratio )

          ratio = 'Unfold / Truth' if config['workflow']=="unfold" else "Reweight / sys. var."
          pltLists[f"Unfoldcomparepseudodata_iter{it}"] = PlotConfig(
                                                                     hcc_data.gen,
                                                                     [hcc_it.gen,hcc_MC.gen],
                                                                     f'pseudodata_truth_unfold_{tag}_iter{it}',
                                                                     ratio)

          ratio = 'Unfold / MC' if config['workflow']=="unfold" else "Reweight / MC"
          pltLists[f"Unfoldcompareeff_iter{it}"]= PlotConfig(
                                                              hcc_MC.eff,
                                                              [hcc_it.eff],
                                                              f'MC_unfoldedeff_{tag}_iter{it}',
                                                              ratio )

          pltLists[f"Unfoldcompareacc_iter{it}"]= PlotConfig( hcc_MC.acc,
                                                                [hcc_it.acc],
                                                                f'MC_unfoldacc_{tag}_iter{it}',
                                                                ratio )

          ratio = 'Unfold / Truth' if config['workflow']=="unfold" else "Reweight / sys. var"
          pltLists[f"Unfoldcomparepseudodataeff_iter{it}"]= PlotConfig(
                                                              hcc_data.eff,
                                                              [hcc_it.eff, hcc_MC.eff],
                                                              f'pseudodatatruth_unfoldedeff_{tag}_iter{it}',
                                                              ratio )

          pltLists[f"Unfoldcomparepseudodataacc_iter{it}"]= PlotConfig( hcc_data.acc,
                                                                [hcc_it.acc,hcc_MC.acc],
                                                                f'pseudodatatruth_unfoldacc_{tag}_iter{it}',
                                                                ratio )
          to_plot = [f'Refoldcompare_iter{it}', f'Unfoldcompare_iter{it}', f'Unfoldcompareeff_iter{it}', f'Unfoldcompareacc_iter{it}',f'Unfoldcomparepseudodata_iter{it}', f"Unfoldcomparepseudodataeff_iter{it}", f"Unfoldcomparepseudodataacc_iter{it}" ]

          print(f"plotting iteration {it}")
          for plt_type in to_plot:
            txt_list = TextListReco if pltLists[plt_type].is_reco else TextListGen
            draw_plot( pltLists[plt_type], config["outputdir"], obs1_name, obs2_name, obs2, txt_list, use_root = False)




