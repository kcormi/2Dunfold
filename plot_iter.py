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
  if use_root:
    plot_flat_hists_root(plt_list.ref.hist, plt_list.hists, plt_list.ref.legend, plt_list.legends, **input_args)
  else:
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
      "merge_bin": True
    }
    for key in kwargs.keys():
      plot_args[key] = kwargs[key]
    if "eff" in plt_list.name:
     plot_args["labelY"] = "Efficiency"
     plot_args["merge_bin"] = False
    elif "acc" in plt_list.name:
     plot_args["labelY"] = "Acceptance"
     plot_args["merge_bin"] = False
    elif "unc_breakdown" in plt_list.name:
     plot_args["labelY"] = "Relative uncertainty" 
     plot_args["merge_bin"] = False
    elif "reco" in plt_list.name:
     plot_args["merge_bin"] = False
    plot_wrapper(plt_list,**plot_args)

    if not("eff" in plt_list.name or "acc" in plt_list.name):
      plot_args["is_logY"] = 1
      plot_wrapper(plt_list,**plot_args)
      

def draw_gof(plt_list, plotdir,**kwargs):
    if plt_list.ref is None: return
    path = f'{plotdir}/{plt_list.name}'
    os.system(f"mkdir -p {plotdir}")
    axis_title = "Iteration"

    plot_args = {
      "title": axis_title,
      "is_logY": 0,
      "do_ratio": 0,
      "output_path": path,
      "text_list": [""],
      "labelY":'Goodness of fit',
      "range_ratio": 0.,
      "use_root": 0
    }
    for key in kwargs.keys():
      plot_args[key] = kwargs[key]
    plot_wrapper(plt_list,**plot_args)
    plot_args["is_logY"] = 1
    plot_wrapper(plt_list,**plot_args)


def get_sel_legend_color(config,config_style,datatype,varname=None):
  if varname is not None:
    sel = {'datatype': "MC",
           'dataset': varname}
    legend = config[f"MCvar_legend"][config[f"MCvar_name"].index(varname)]
    color = config_style['mpl']["MCvar_color"][config[f"MCvar_name"].index(varname)]
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
      result_settings.gen_reweight = (config['workflow']=='genreweight')
      result_settings.sys_reweight = (config['workflow']=='sysreweight')

      chi2_collection_list = []
      chi2_collection_unfold_MC_list = []
      chi2_collection_MC_data_list = []
      ks_collection_list = []
      for icolor,obs in enumerate(args.obs):
      
        (obs1_name, obs1_binning), (obs2_name, obs2_binning) = parse_obs_args( config, obs )

        obs1 = ObsConfig.from_yaml( config["varunfold"], [obs1_name], binning=obs1_binning )
        obs2 = ObsConfig.from_yaml( config["varunfold"], [obs2_name], binning=obs2_binning )
        TextListReco = [f'{obs1.reco.edges[i]} #leq {obs1.reco.shortname} < {obs1.reco.edges[i+1]}' for i in range(obs1.reco.nbins)]
        TextListGen = [f'{obs1.gen.edges[i]} #leq {obs1.gen.shortname} < {obs1.gen.edges[i+1]}' for i in range(obs1.gen.nbins)]+(["background"] if config["addbkg"] else [])

        base_sel = {'obs1': obs[0],
                    'obs2': obs[1],
                    'method': method
                   }
        if "pseudodata" in config.keys() and config["pseudodata"]:
          data_sel, data_legend, data_color = get_sel_legend_color(config,config_style,"pseudodata")
          chi2_sel = {'datatype': f"chi2_{config['workflow']}_pseudodata",
                      'dataset':f"{config['workflow']}_{config['pseudodata_name']}"}
          ks_sel = {'datatype': f"ks_{config['workflow']}_pseudodata",
                      'dataset':f"{config['workflow']}_{config['pseudodata_name']}"}
          chi2_MC_data_sel = {'datatype': f"chi2_MC_pseudodata",
                      'dataset':f"{config['MC_name']}_{config['pseudodata_name']}"}
        else:
          data_sel, data_legend, data_color = get_sel_legend_color(config,config_style,"data")
          chi2_sel = {'datatype': f"chi2_{config['workflow']}_data",
                      'dataset':f"{config['workflow']}_{config['data_name']}"}
          ks_sel = {'datatype': f"ks_{config['workflow']}_data",
                      'dataset':f"{config['workflow']}_{config['data_name']}"}
          chi2_MC_data_sel = {'datatype': f"chi2_MC_data",
                      'dataset':f"{config['MC_name']}_{config['data_name']}"}
          print("chi2_unfold_MC_sel",chi2_MC_data_sel)
        chi2_unfold_MC_sel = {'datatype': f"chi2_{config['workflow']}_MC",
                    'dataset':f"{config['workflow']}_{config['MC_name']}"}
        print("chi2_unfold_MC_sel",chi2_unfold_MC_sel)
        MC_sel, MC_legend, MC_color = get_sel_legend_color(config,config_style,"MC")
        MCvar_sel = {}
        MCvar_legend = {}
        MCvar_color = {}
        if "MCvar_name" in config.keys():
          for MCvar in config["MCvar_name"]:
            MCvar_sel[MCvar], MCvar_legend[MCvar], MCvar_color[MCvar] = get_sel_legend_color(config,config_style,"MC",varname=MCvar)


        records_data = get_df_entries(df, **base_sel|data_sel)
        records_MC = get_df_entries(df, **base_sel|MC_sel)
        records_chi2 = get_df_entries(df, **base_sel|chi2_sel)
        records_ks = get_df_entries(df, **base_sel|ks_sel)

        print(base_sel)
        records_chi2_unfold_MC = get_df_entries(df, **base_sel|chi2_unfold_MC_sel)
        records_chi2_MC_data = get_df_entries(df, **base_sel|chi2_MC_data_sel)
        print("records_chi2_unfold_MC",records_chi2_unfold_MC)
        print("records_chi2_MC_data",records_chi2_MC_data)

        hcc_data = get_histconfigcollection(records_data,data_color,['boldertriangle','boldercross','boldertriangle','boldertriangle',None],data_legend,norm = (False if ("normalize_XS" in config.keys() and config["normalize_XS"]) else True))

        hcc_MC = get_histconfigcollection(records_MC,MC_color,['fillederror','fillederror','cross','cross',None],MC_legend,norm = (False if ("normalize_XS" in config.keys() and config["normalize_XS"]) else True))
        hcc_MCvar_list = []
        if "MCvar_name" in config.keys():
          for MCvar in config["MCvar_name"]:
            hcc_MCvar_list.append(get_histconfigcollection(get_df_entries(df, **base_sel|MCvar_sel[MCvar]),MCvar_color[MCvar],['fillederror','fillederror','cross','cross',None],MCvar_legend[MCvar],norm = (False if ("normalize_XS" in config.keys() and config["normalize_XS"]) else True)))


        chi2_collection,_ = get_histconfig_gof(records_chi2,default_colors[icolor],(f"{obs[0]}:{obs[1]}").replace("spherocity","sphericity"),["gen_both_error","gen_target_error","reco_both_error","reco_target_error"])
        chi2_collection_list.append(chi2_collection)

        chi2_collection_unfold_MC,_ = get_histconfig_gof(records_chi2_unfold_MC,default_colors[icolor],obs[1].split(":")[0].replace("spherocity","sphericity"),["gen_both_error","gen_target_error","reco_both_error","reco_target_error"])
        chi2_collection_MC_data,_ = get_histconfig_gof(records_chi2_MC_data,default_colors[icolor],(f"{obs[0]}:{obs[1]}").replace("spherocity","sphericity"),["gen_both_error","gen_target_error","reco_both_error","reco_target_error"])
        print("records_chi2_MC_data",records_chi2_MC_data)
        print("chi2_collection_MC_data",chi2_collection_MC_data)
        chi2_collection_unfold_MC_list.append(chi2_collection_unfold_MC)
        chi2_collection_MC_data_list.append(chi2_collection_MC_data)


        _,ks_collection = get_histconfig_gof(records_ks,default_colors[icolor],(f"{obs[0]}:{obs[1]}").replace("spherocity","sphericity"),["gen_distance","reco_distance","gen_p","reco_p"])
        ks_collection_list.append(ks_collection)

        iters = np.sort(np.unique(df[~np.isnan(df['iter'])]['iter']))
       
        for it in iters:
          iter_sel = {'datatype': config['workflow']+args.unfold_tag,
                      'dataset': config['MC_name'],
                      'iter': it}
          records_it = get_df_entries(df, **base_sel|iter_sel)
          it_color = result_settings.color
          it_color_unfold = result_settings.color_unfold
          it_legend = result_settings.legend_refold+args.unfold_tag
          it_legend_unfold = result_settings.legend_unfold+args.unfold_tag
          tag = result_settings.tag
          hcc_it = get_histconfigcollection(records_it,[it_color_unfold,it_color,it_color_unfold,it_color_unfold,None],["marker","filled","cross","cross",None],[it_legend_unfold,it_legend,it_legend_unfold+" Eff.",it_legend_unfold+" Acc.",it_legend_unfold+" Mig."],norm = (False if ("normalize_XS" in config.keys() and config["normalize_XS"]) else True))

          hcc_sys_list = []
          if method in config["evaluate-relunc"]:
            hc_relunc_list = []
            ha_relunc_dic = {}
            color_unc = ['k','c','m','y','steelblue','r','b','indianred','saddlebrown']
            for isys,sys in enumerate(config["reluncname"]):
              if sys != "totalsys" and sys != "MCstat":
                sys_sel = {"datatype": config['workflow'],
                           "dataset": sys,
                           'iter': it
                           }
                records_sys = get_df_entries(df, **base_sel|sys_sel)
                print(base_sel|sys_sel)
                it_legend_sys = config["relunclegend"][isys]
                hcc_sys_list.append(get_histconfigcollection(records_sys,[color_unc[isys],color_unc[isys],color_unc[isys],color_unc[isys],None],["marker","fillederror","cross","cross",None],[it_legend_sys,it_legend_sys,it_legend_sys+" Eff.",it_legend_sys+" Acc.",it_legend_sys+" Mig."],norm = (False if ("normalize_XS" in config.keys() and config["normalize_XS"]) else True)))
                print(records_sys)
                histarray_sys = HistArray.from_df_entry(records_sys.iloc[0])
                histarray_sys.divide_by_bin_width()
                ha_relunc_dic[sys] = abs(hcc_sys_list[-1].gen.hist - hcc_it.gen.hist)/hcc_it.gen.hist
            if "MCstat" in config["reluncname"]:
              ha_relunc_dic["MCstat"] = hcc_it.gen.hist.relative_error()
            if "totalsys" in config["reluncname"]:
              ha_relunc_dic["totalsys"] = sum([ha_relunc_dic[sys_decompose]**2 for sys_decompose in ha_relunc_dic.keys()])**(1./2.)
            for isys,sys in enumerate(config["reluncname"]): 
              hc_relunc_list.append(HistConfig(ha_relunc_dic[sys],0,
                                        color_unc[isys], "line",config["relunclegend"][isys]))

            

          if method in config['evaluate-sys-bootstrap']:
            iter_sysbs_sel = {'datatype': config['workflow']+"_sys_fixXS",
                      'dataset': config['MC_name'],
                      'iter': it}
            records_it_sysbs = get_df_entries(df, **base_sel|iter_sysbs_sel)
            it_sysbs_color = result_settings.color
            it_sysbs_color_unfold = result_settings.color_unfold
            it_sysbs_legend = result_settings.legend_refold+" sys."
            it_sysbs_legend_unfold = result_settings.legend_unfold+" sys."
            hcc_it_sysbs = get_histconfigcollection(records_it_sysbs,[it_sysbs_color_unfold,it_sysbs_color,it_sysbs_color_unfold,it_sysbs_color_unfold,None],["hatch","hatch","hatch","hatch",None],[it_sysbs_legend_unfold,it_sysbs_legend,it_sysbs_legend_unfold+" Eff.",it_sysbs_legend_unfold+" Acc.",it_sysbs_legend_unfold+" Mig."],norm = (False if ("normalize_XS" in config.keys() and config["normalize_XS"]) else True))
            if method in config["evaluate-relunc"]:
              hc_relunc_list.append(HistConfig(hcc_it_sysbs.gen.hist.relative_error(),0,
                                        'k', "dashedline","Total sys. unc. (pseudo-experiment)"))
          else:
            hcc_it_sysbs = None

          if method in config['evaluate-stat-bootstrap']:
            iter_statbs_sel = {'datatype': config['workflow']+"_stat_fixXS",
                      'dataset': config['MC_name'],
                      'iter': it}
            records_it_statbs = get_df_entries(df, **base_sel|iter_statbs_sel)
            it_statbs_color = result_settings.color
            it_statbs_color_unfold = result_settings.color_unfold
            it_statbs_legend = result_settings.legend_refold+" stat."
            it_statbs_legend_unfold = result_settings.legend_unfold+" stat."
            hcc_it_statbs = get_histconfigcollection(records_it_statbs,[it_statbs_color_unfold,it_statbs_color,it_statbs_color_unfold,it_statbs_color_unfold,None],["righthatch","righthatch","righthatch","righthatch",None],[it_statbs_legend_unfold,it_statbs_legend,it_statbs_legend_unfold+" Eff.",it_statbs_legend_unfold+" Acc.",it_statbs_legend_unfold+" Mig."],norm = (False if ("normalize_XS" in config.keys() and config["normalize_XS"]) else True))
            if method in config["evaluate-relunc"]:
              hc_relunc_list.append(HistConfig(hcc_it_statbs.gen.hist.relative_error(),0,
                                        'orange', "dashedline","Stat. unc. (pseudo-experiment)"))
          else:
            hcc_it_statbs = None

          if method in config['evaluate-sys-bootstrap'] and method in config['evaluate-stat-bootstrap']:
            iter_fullbs_sel = {'datatype': config['workflow']+"_fullbsunc_fixXS",
                      'dataset': config['MC_name'],
                      'iter': it}
            records_it_fullbs = get_df_entries(df, **base_sel|iter_fullbs_sel)
            it_fullbs_color = result_settings.color
            it_fullbs_color_unfold = "black"
            it_fullbs_legend = result_settings.legend_refold
            it_fullbs_legend_unfold = result_settings.legend_unfold
            hcc_it_fullbs = get_histconfigcollection(records_it_fullbs,[it_fullbs_color_unfold,it_fullbs_color,it_fullbs_color_unfold,it_fullbs_color_unfold,None],["marker","marker","marker","marker",None],[it_fullbs_legend_unfold,it_fullbs_legend,it_fullbs_legend_unfold+" Eff.",it_fullbs_legend_unfold+" Acc.",it_fullbs_legend_unfold+" Mig."],records_it_statbs,norm = (False if ("normalize_XS" in config.keys() and config["normalize_XS"]) else True))
          else:
            hcc_it_fullbs = None



          pltLists = {}
          ratio = 'Refold / data' if config['workflow']=="unfold" else "Reweight / sys. var."
          pltLists[f"Refoldcompare_iter{it}"] = PlotConfig(
                                                           hcc_data.reco,
                                                           [hcc_it.reco,hcc_it_sysbs.reco if hcc_it_sysbs else None,hcc_it_statbs.reco if hcc_it_statbs else None,hcc_MC.reco]+[hcc_sys.reco for hcc_sys in hcc_sys_list],
                                                           f'data_refold_{tag}_iter{it}',
                                                           ratio)

          ratio = 'Unfold / MC' if config['workflow']=="unfold" else "Reweight / MC"
          pltLists[f"Unfoldcompare_iter{it}"] = PlotConfig(
                                                           hcc_MC.gen,
                                                           [hcc_it.gen,hcc_it_sysbs.gen if hcc_it_sysbs else None,hcc_it_statbs.gen if hcc_it_statbs else None]+[hcc_sys.gen for hcc_sys in hcc_sys_list],
                                                           f'MC_unfold_{tag}_iter{it}',
                                                           ratio )

          if method in config['evaluate-sys-bootstrap'] and method in config['evaluate-stat-bootstrap']:
            ratio = 'MC / Data'
            pltLists[f"Unfoldbscompare_iter{it}"] = PlotConfig(
                                                             hcc_it_fullbs.gen,
                                                             [hcc_MC.gen]+[hcc_MCvar.gen for hcc_MCvar in hcc_MCvar_list],
                                                             f'MC_unfold_bsunc_{tag}_iter{it}',
                                                             ratio)
          
          ratio = 'Unfold / Truth' if config['workflow']=="unfold" else "Reweight / sys. var."
          pltLists[f"Unfoldcomparepseudodata_iter{it}"] = PlotConfig(
                                                                     hcc_data.gen,
                                                                     [hcc_it.gen,hcc_it_sysbs.gen if hcc_it_sysbs else None,hcc_it_statbs.gen if hcc_it_statbs else None,hcc_MC.gen]+[hcc_sys.gen for hcc_sys in hcc_sys_list]+[hcc_MCvar.gen for hcc_MCvar in hcc_MCvar_list],
                                                                     f'pseudodata_truth_unfold_{tag}_iter{it}',
                                                                     ratio)


            

          ratio = 'Unfold / MC' if config['workflow']=="unfold" else "Reweight / MC"
          pltLists[f"Unfoldcompareeff_iter{it}"]= PlotConfig(
                                                              hcc_MC.eff,
                                                              [hcc_it.eff,hcc_it_sysbs.eff if hcc_it_sysbs else None,hcc_it_statbs.eff if hcc_it_statbs else None]+[hcc_sys.eff for hcc_sys in hcc_sys_list],
                                                              f'MC_unfoldedeff_{tag}_iter{it}',
                                                              ratio )

          pltLists[f"Unfoldcompareacc_iter{it}"]= PlotConfig( hcc_MC.acc,
                                                                [hcc_it.acc,hcc_it_sysbs.acc if hcc_it_sysbs else None,hcc_it_statbs.acc if hcc_it_statbs else None]+[hcc_sys.acc for hcc_sys in hcc_sys_list],
                                                                f'MC_unfoldacc_{tag}_iter{it}',
                                                                ratio )

          ratio = 'Unfold / Truth' if config['workflow']=="unfold" else "Reweight / sys. var"
          pltLists[f"Unfoldcomparepseudodataeff_iter{it}"]= PlotConfig(
                                                              hcc_data.eff,
                                                              [hcc_it.eff,hcc_it_sysbs.eff if hcc_it_sysbs else None,hcc_it_statbs.eff if hcc_it_statbs else None, hcc_MC.eff]+[hcc_sys.eff for hcc_sys in hcc_sys_list],
                                                              f'pseudodatatruth_unfoldedeff_{tag}_iter{it}',
                                                              ratio )

          pltLists[f"Unfoldcomparepseudodataacc_iter{it}"]= PlotConfig( hcc_data.acc,
                                                                [hcc_it.acc,hcc_it_sysbs.acc if hcc_it_sysbs else None,hcc_it_statbs.acc if hcc_it_statbs else None,hcc_MC.acc]+[hcc_sys.acc for hcc_sys in hcc_sys_list],
                                                                f'pseudodatatruth_unfoldacc_{tag}_iter{it}',
                                                                ratio )
          if method in config["evaluate-relunc"]:
            legend = "Unc./Total"
            if method in config['evaluate-stat-bootstrap'] or method in config['evaluate-sys-bootstrap']: legend = "Unc./Total sys.(quadrature sum)"
            pltLists[f"RelativeUncertainty_iter{it}"] = PlotConfig(hc_relunc_list[0],
                                                                  hc_relunc_list[1:] if len(hc_relunc_list)>1 else [],
                                                                  f'unc_breakdown_{tag}_iter{it}',legend)

          
          to_plot = [f'Refoldcompare_iter{it}', f'Unfoldcompare_iter{it}', f'Unfoldcompareeff_iter{it}', f'Unfoldcompareacc_iter{it}',f'Unfoldcomparepseudodata_iter{it}', f"Unfoldcomparepseudodataeff_iter{it}", f"Unfoldcomparepseudodataacc_iter{it}"]+([f"RelativeUncertainty_iter{it}" ] if method in config["evaluate-relunc"] else []) + ([f"Unfoldbscompare_iter{it}"] if method in config['evaluate-sys-bootstrap'] and method in config['evaluate-stat-bootstrap'] else [])
          print(f"plotting iteration {it}")
          for plt_type in to_plot:
            txt_list = TextListReco if pltLists[plt_type].is_reco else TextListGen
            print(plt_type,"drawing")
            draw_plot( pltLists[plt_type], config["outputdir"], obs1_name, obs2_name, obs2, txt_list, use_root = False, normalize_XS = True if ("normalize_XS" in config.keys() and config["normalize_XS"]) else False)
            print("draw finished")
        
      chi2_plot = ["gen_target_error","gen_target_error_rel","gen_both_error","gen_both_error_rel","reco_target_error","reco_target_error_rel","reco_both_error","reco_both_error_rel"]
      #ks_plot = ["gen_distance","gen_distance_rel","reco_distance","reco_distance_rel","gen_p","gen_p_rel","reco_p","reco_p_rel"]
      ks_plot = ["gen_distance","gen_distance_rel","reco_distance","reco_distance_rel"]
      for (gof,gof_plot,gof_collection_list) in zip(["chi2","ks"],[chi2_plot,ks_plot],[chi2_collection_list,ks_collection_list]):
        for name in gof_plot:
          gof_list = [getattr(gof_coll,name) for gof_coll in gof_collection_list]
          gof_list = [histconfig for histconfig in gof_list if histconfig is not None]
          if len(gof_list)>0:
            gof_mean = copy.deepcopy(gof_list[0])
            gof_mean.color='k'
            gof_mean.legend = 'mean'
            gof_mean.hist = sum([gof_obs.hist for gof_obs in gof_list])
            gof_mean.hist = gof_mean.hist/len(gof_list)
            gof_list.insert(0,gof_mean)
            plotconfig = PlotConfig( gof_list[0],
                                     gof_list[1:] if len(gof_list)>1 else [],
                                     f"{gof}_{name}_{method}",
                                     "" )
            draw_gof( plotconfig, config["outputdir"] )
            print("draw finished")
        
        
      chi2_plot_bottomline = ["both_error_bottomline"]
      chi2_plot_unfold_MC = ["gen_both_error"]
      chi2_plot_MC_data = ["reco_both_error"]
      chi2_plot_refold_MC = ["reco_both_error"]
      print("chi2_collection_unfold_MC_list",chi2_collection_unfold_MC_list)
      print("chi2_collection_MC_data_list",chi2_collection_MC_data_list)
      for (name_bottomline,name_unfold_MC,name_MC_data) in zip(chi2_plot_bottomline,chi2_plot_unfold_MC,chi2_plot_MC_data):
        gof_unfold_MC_list = [getattr(gof_coll,name_unfold_MC) for gof_coll in chi2_collection_unfold_MC_list]
        gof_unfold_MC_list = [histconfig for histconfig in gof_unfold_MC_list if histconfig is not None]
        gof_MC_data_list = [getattr(gof_coll,name_MC_data) for gof_coll in chi2_collection_MC_data_list]
        gof_MC_data_list = [histconfig for histconfig in gof_MC_data_list if histconfig is not None]
        if len(gof_unfold_MC_list)>0 and len(gof_MC_data_list)>0:
          print(name_bottomline,"drawing")
          gof_bottomline_list=[]
          for(gof_unfold_MC,gof_MC_data) in zip(gof_unfold_MC_list,gof_MC_data_list):
              gof_bottomline = copy.deepcopy(gof_unfold_MC)
              #print("gof_MC_data.hist.nested_value",gof_MC_data.hist.nested_value)
              gof_bottomline.hist = gof_unfold_MC.hist/gof_MC_data.hist.nested_value[0][0]
              gof_bottomline_list.append(gof_bottomline)
          plotconfig = PlotConfig( gof_bottomline_list[0],
                                       gof_bottomline_list[1:] if len(gof_bottomline_list)>1 else [],
                                       f"chi2_{name_bottomline}_{method}",
                                       "" )
          draw_gof( plotconfig, config["outputdir"],labelY=r'$\frac{\chi^2(\mathrm{unfold & gen MC)}}{\chi^2(\mathrm{data & smeared MC})}$')
        
      for (name_bottomline,name_unfold_MC,name_refold_MC) in zip(chi2_plot_bottomline,chi2_plot_unfold_MC,chi2_plot_refold_MC):
        gof_unfold_MC_list = [getattr(gof_coll,name_unfold_MC) for gof_coll in chi2_collection_unfold_MC_list]
        gof_unfold_MC_list = [histconfig for histconfig in gof_unfold_MC_list if histconfig is not None]
        gof_refold_MC_list = [getattr(gof_coll,name_refold_MC) for gof_coll in chi2_collection_unfold_MC_list]
        gof_refold_MC_list = [histconfig for histconfig in gof_refold_MC_list if histconfig is not None]
        if len(gof_unfold_MC_list)>0 and len(gof_refold_MC_list)>0:
          print(name_bottomline,"drawing")
          gof_bottomline_list=[]
          for(gof_unfold_MC,gof_refold_MC) in zip(gof_unfold_MC_list,gof_refold_MC_list):
            gof_bottomline = copy.deepcopy(gof_unfold_MC)
            #print("gof_MC_data.hist.nested_value",gof_MC_data.hist.nested_value)
            gof_bottomline.hist = gof_unfold_MC.hist/gof_refold_MC.hist
            gof_bottomline_list.append(gof_bottomline)
          plotconfig = PlotConfig( gof_bottomline_list[0],
                                     gof_bottomline_list[1:] if len(gof_bottomline_list)>1 else [],
                                     f"chi2_{name_bottomline}_refold_{method}",
                                     "" )
          draw_gof( plotconfig, config["outputdir"],labelY=r'$\frac{\chi^2(\mathrm{unfold & gen MC)}}{\chi^2(\mathrm{refold & smeared MC})}$' )
        
