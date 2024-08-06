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
from ast import literal_eval
import copy
@dataclass
class HistConfigCollection:
  gen: HistConfig = None
  reco: HistConfig = None
  eff: HistConfig = None
  acc: HistConfig = None
  mig: HistConfig = None


@dataclass
class Chi2Collection:
  gen_both_error: HistConfig = None
  gen_both_error_rel: HistConfig = None
  gen_target_error: HistConfig = None
  gen_target_error_rel: HistConfig = None
  reco_both_error: HistConfig = None
  reco_both_error_rel: HistConfig = None
  reco_target_error: HistConfig = None
  reco_target_error_rel: HistConfig = None

@dataclass
class KSCollection:
  gen_distance: HistConfig = None
  gen_distance_rel: HistConfig = None
  reco_distance: HistConfig = None
  reco_distance_rel: HistConfig = None
  gen_p: HistConfig = None
  gen_p_rel: HistConfig = None
  reco_p: HistConfig = None 
  reco_p_rel: HistConfig = None



def get_df_entries(df, **kwargs_selections):
  for ikey,key in enumerate(kwargs_selections.keys()):
    if ikey == 0:
      selection = (df[key]==kwargs_selections[key])
    else:
      selection = selection & (df[key]==kwargs_selections[key])
    
  return df[selection]

def get_histconfigcollection(df,color_list,style_list,legend_list,df_stat=None,norm=False):
  gen_entry = df[(df['histtype']=='inclusive')&df['dim1_isgen']&df['dim2_isgen']]
  reco_entry = df[(df['histtype']=='inclusive')&~df['dim1_isgen']&~df['dim2_isgen']] 
  #Plan: process the results from bootstraps -> add 'histtype' options 'inclusive_sysunc', 'inclusive_statunc','inclusive_totalunc'
  eff_entry = df[(df['histtype']=='efficiency')&df['dim1_isgen']&df['dim2_isgen']]
  acc_entry = df[(df['histtype']=='acceptance')&~df['dim1_isgen']&~df['dim2_isgen']]
  mig_entry = df[(df['histtype']=='migration')&df['dim1_isgen']&~df['dim2_isgen']] 

  entries = [gen_entry,reco_entry,eff_entry,acc_entry,mig_entry]

  if not isinstance(color_list,list):
    color_list = [color_list]*4+[None]

  if not isinstance(style_list,list):
    style_list = [style_list]*4+[None]

  if not isinstance(legend_list,list):
    legend_list = [legend_list,legend_list,legend_list+" Eff.",legend_list+" Acc.",legend_list+" Mig."]

  entries_statonly = [0]*5
  if df_stat is not None:
    entries_statonly = [df_stat[(df_stat['histtype']=='inclusive')&df_stat['dim1_isgen']&df_stat['dim2_isgen']],
                          df_stat[(df_stat['histtype']=='inclusive')&~df_stat['dim1_isgen']&~df_stat['dim2_isgen']],
                          df_stat[(df_stat['histtype']=='efficiency')&df_stat['dim1_isgen']&df_stat['dim2_isgen']],
                          df_stat[(df_stat['histtype']=='acceptance')&~df_stat['dim1_isgen']&~df_stat['dim2_isgen']],
                          df_stat[(df_stat['histtype']=='migration')&df_stat['dim1_isgen']&~df_stat['dim2_isgen']]]

  norm=False
  list_HistConfig = []
  for i,(entry,entry_statonly) in enumerate(zip(entries,entries_statonly)):
    if len(entry)>0:
      histarray = HistArray.from_df_entry(entry.iloc[0])
      if i<2:
        if norm:
          histarray.normalize()
        histarray.divide_by_bin_width()
      if df_stat is not None and len(entry_statonly)>0:
        histarray_statonly = HistArray.from_df_entry(entry_statonly.iloc[0])
        if i<2:
          if norm:
            histarray_statonly.normalize()
          histarray_statonly.divide_by_bin_width()
      else:
        histarray_statonly = 0
      list_HistConfig.append(HistConfig(histarray,
                                        histarray_statonly,
                                        color_list[i], style_list[i],legend_list[i]))
    else:
      list_HistConfig.append(None)

  return HistConfigCollection(*list_HistConfig)

def get_histconfig_gof(df,color,legend,histtype_list):
  #gen_both_error_entry = df[(df['histtype']=='gen_both_error')]
  #gen_target_error_entry = df[(df['histtype']=='gen_target_error')]
  #reco_both_error_entry = df[(df['histtype']=='reco_both_error')]
  #reco_target_error_entry = df[(df['histtype']=='reco_target_error')]

  list_gof = []
  #for i,entry in enumerate([gen_both_error_entry,gen_target_error_entry,reco_both_error_entry,reco_target_error_entry]):
  for i,histtype in enumerate(histtype_list):
    entry = df[(df['histtype']==histtype)]
    if len(entry)>0:
      histarray = HistArray.from_df_entry(entry.iloc[0])
      list_gof.append(HistConfig(histarray,0,color,"hist",legend))
      list_gof.append(HistConfig(histarray/histarray.nested_value[0][0],0,color,"hist",legend))
    else:
      list_gof.append(None)
      list_gof.append(None)

  return Chi2Collection(*list_gof),KSCollection(*list_gof)


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
  plot_flat_hists_mpl(plt_list.ref.hist, plt_list.hists, plt_list.ref.legend, plt_list.legends, **input_args)

def draw_plot(plt_list, plotdir, var1_nm, var2_nm, obs2, txt_list,use_root=True,normalize_XS=True,**kwargs):
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
      "labelY":r'$\sigma$/d'+latex_root_to_mpl(axis_title) if normalize_XS else "1/N dN/d"+r'{}'.format(latex_root_to_mpl(axis_title).replace(']','^{-1}]')),
      "range_ratio": 0.1 if "eff" in plt_list.name or "acc" in plt_list.name else 0.3,
      "use_root": use_root,
      "merge_bin": False
    }
    for key in kwargs.keys():
      plot_args[key] = kwargs[key]
    if "reco" in plt_list.name:
     plot_args["merge_bin"] = False
    plot_wrapper(plt_list,**plot_args)

    if not("eff" in plt_list.name or "acc" in plt_list.name):
      plot_args["is_logY"] = 1
      plot_wrapper(plt_list,**plot_args)


def get_sel_legend_color(config,config_style,datatype,varname=None):
  if varname is not None:
    sel = {'datatype': "MC",
           'dataset': varname}
  else:
    sel = {'datatype': datatype,
           'dataset': config[f"{datatype}_name"]}
    legend = config[f"{datatype}_legend"]
    color = config_style['mpl'][f"{datatype}_color"]
  return sel,legend,color

if __name__=="__main__":

    parser = ArgumentParser()
    parser.add_argument('--config',default="config/plot_1d_v7_MCEPOS_unfoldCP1.json",help="The configration file including the unfolding setup")
    parser.add_argument('--config-style', type=str, default = 'config/results_style.yml')
    parser.add_argument('--obs', type=obs_pair,nargs = '+', help="Pair of observables separated by a comma")
    parser.add_argument('--unfold-tag', type=str, default = '',help="To check the unfolding results of bootstraps,  set it into e.g. _sysbs_1")
    args = parser.parse_args()

    if len(args.obs) < 11:
      default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    else:
      cmap = plt.get_cmap('jet')
      default_colors = cmap(np.linspace(0, 1.0, len(args.obs)))

    with open(args.config, 'r') as configjson:
        config = json.load(configjson)

    style_file = args.config_style
    with open(style_file, 'r') as style_json:
        config_style = yaml.safe_load(style_json)

    df = pd.read_csv(config["input"])

    for method in config["methods"]:
      result_settings = ResultPlotSettings.from_yaml(style_file, ['mpl', method])

      for icolor,obs in enumerate(args.obs):
        (obs1_name, obs1_binning), (obs2_name, obs2_binning) = parse_obs_args( config, obs )


        base_sel = {'obs1': obs[0],
                    'obs2': obs[1],
                    'method': method,
                    'histtype': 'inclusive'
                   }
        data_sel, data_legend, data_color = get_sel_legend_color(config,config_style,"data")
        MC_sel, MC_legend, MC_color = get_sel_legend_color(config,config_style,"MC")
        if "MC_color" in config:
            MC_color = config["MC_color"]

        print(df)
        print(df[ df['datatype'] == 'data' ])
        print(base_sel, data_sel)
        records_data = get_df_entries(df, **base_sel|data_sel)
        print(f'data records: {records_data}')
        print(f'data records: {records_data.iloc[0].bin_edges_dim1}')
        print(f'data records: {records_data.iloc[0].bin_edges_dim2}')
        #print(f'data records: {records_data.iloc[0].bin_values}')
        records_MC = get_df_entries(df, **base_sel|MC_sel)

        outer_edges = literal_eval(records_data.iloc[0].bin_edges_dim1)
        obs1 = ObsConfig.from_yaml( config["varunfold"], [obs1_name], binning=obs1_binning )
        obs2 = ObsConfig.from_yaml( config["varunfold"], [obs2_name], binning=obs2_binning )
        TextListReco = [f'{round(outer_edges[i]):d} #leq {obs1.reco.shortname} < {round(outer_edges[i+1]):d}' for i, _ in enumerate(outer_edges[:-1])]
        TextListGen = [f'{obs1.gen.edges[i]} #leq {obs1.gen.shortname} < {obs1.gen.edges[i+1]}' for i in range(obs1.gen.nbins)]+(["background"] if config["addbkg"] else [])

        print(f'base selection: {base_sel}')
        hcc_data = get_histconfigcollection(records_data,data_color,['boldertriangle','boldercross','boldertriangle','boldertriangle',None],data_legend,norm = (False if ("normalize_XS" in config.keys() and config["normalize_XS"]) else True))
        #print(f'{hcc_data} values, edges: {hcc_data.reco.hist.nested_value}, {hcc_data.reco.hist.nested_bins}')

        hcc_MC = get_histconfigcollection(records_MC,MC_color,['marker','fillederror','cross','cross',None],MC_legend,norm = (False if ("normalize_XS" in config.keys() and config["normalize_XS"]) else True))

        for it in ['dummy']:
          iter_sel = {'datatype': config['workflow']+args.unfold_tag,
                      'dataset': config['MC_name']}
          records_it = get_df_entries(df, **base_sel|iter_sel)
          it_color = result_settings.color
          it_color_unfold = result_settings.color_unfold
          it_legend = result_settings.legend_refold+args.unfold_tag
          it_legend_unfold = result_settings.legend_unfold+args.unfold_tag
          tag = result_settings.tag

          hcc_sys_list = []
          if method in config["evaluate-relunc"]:
            #color_unc = ['orchid','darkviolet','darkmagenta','c','m','y','steelblue','r','b','indianred','saddlebrown']
            color_unc = ['darkviolet','darkmagenta','c','m','y','steelblue','r','b','indianred','saddlebrown']
            for isys,sys in enumerate(config["reluncname"]):
                print(f'Running {sys}')
                sys_sel = {"datatype": config['workflow'],
                           "dataset": sys,
                            "dim1_isgen": False,
                           }
                records_sys = get_df_entries(df, **base_sel|sys_sel)
                #print(records_sys)
                it_legend_sys = config["relunclegend"][isys]
                color = color_unc[isys] if not "relunccol" in config else config["relunccol"][isys]
                style = "fillederror" if not sys.startswith('instanton') and not 'fitDiag' in sys else "hist"
                if "reluncstyle" in config:
                    style = config["reluncstyle"][isys]
                hcc_sys_list.append(get_histconfigcollection(records_sys,[color,color,color,color,None],["marker",style,"cross","cross",None],[it_legend_sys,it_legend_sys],norm = (False if ("normalize_XS" in config.keys() and config["normalize_XS"]) else True)))
                histarray_sys = HistArray.from_df_entry(records_sys.iloc[0])
                histarray_sys.divide_by_bin_width()

          hcc_it_sysbs = None
          hcc_it_statbs = None
          hcc_it_fullbs = None

          #print( hcc_MC.reco.hist.nested_value, hcc_MC.reco.hist.nested_bins)
          pltLists = {}
          ratio = 'pred. / data' 
          pltLists[f"Refoldcompare_"] = PlotConfig(
                                                           hcc_data.reco,
                                                           [None if True else hcc_it.reco, None, None,hcc_MC.reco]+[hcc_sys.reco for hcc_sys in hcc_sys_list],
                                                           f'data_refold_{tag}',
                                                           ratio)

          #ratio = 'Unfold / MC' if config['workflow']=="unfold" else "Reweight / MC"
          #pltLists[f"Unfoldcompare_iter{it}"] = PlotConfig(
          #                                                 hcc_MC.gen,
          #                                                 [hcc_it.gen,None,None]+[hcc_sys.gen for hcc_sys in hcc_sys_list],
          #                                                 f'MC_unfold_{tag}_iter{it}',
          #                                                 ratio )


          to_plot = [f'Refoldcompare_']#, f'Unfoldcompare_iter{it}']#,f'Unfoldcomparepseudodata_iter{it}']
          plot_kwargs= {}
          if "ratio_limits" in config:
                plot_kwargs.update({"ratio_limits":config["ratio_limits"]})
          if "ratio_type" in config:
                plot_kwargs.update({"ratio_type":config["ratio_type"]})
          if "limY" in config:
                plot_kwargs.update({"limY": config["limY"]})
          for plt_type in to_plot:
            txt_list = TextListReco if pltLists[plt_type].is_reco else TextListGen
            print(plt_type,"drawing")
            print([(h.nested_value, h.nested_bins) for h in pltLists[plt_type].hists])
            draw_plot( pltLists[plt_type], config["outputdir"], obs1_name, obs2_name, obs2, txt_list, use_root = False, normalize_XS = True if ("normalize_XS" in config.keys() and config["normalize_XS"]) else False, **plot_kwargs)
            print(f"draw finished: plots written two {config['outputdir']}")
