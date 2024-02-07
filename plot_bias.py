import numpy as np
import json
import yaml
from argparse import ArgumentParser
import os
import ROOT
import re
from unfold_utils import *
from arg_parsing import *
from configs import ObsConfig
import pandas as pd
from GOF.binned import *
import uproot
import awkward as ak
import subprocess
import ast
from multiprocessing.pool import Pool
from itertools import combinations_with_replacement
import copy
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")
import seaborn as sns

ONE_SIGMA = 0.682
ONE_SIGMA_PERCS = [ 100*(1 - ONE_SIGMA)/2, 100*((1 + ONE_SIGMA)/2) ]


def get_bs_entries(df,label_row):
  for ikey, key in enumerate(list(label_row.keys())):
    if ikey == 0:
      sel = (df[key]==label_row[key])
    else:
      sel = sel & (df[key]==label_row[key])
  return df[sel]

if __name__=="__main__":

    parser = ArgumentParser()
    parser.add_argument('-i','--input', default="results_finebin_v7_MCA3P_unfoldCUEP8M1_maxweight10_minweightp1_coverage_v3_1/merge_bs_dataunc.csv", help="The .csv file storing the bias and coverage")
    parser.add_argument('-o','--output', default="results_finebin_v7_MCA3P_unfoldCUEP8M1_maxweight10_minweightp1_coverage_v3_1/plots_multifold_bias_dataunc/", help="The directory for bias plots")
    args = parser.parse_args()
    os.system("mkdir -p "+args.output)

    labels_obs = ['obs1','obs2','iter']
    df = pd.read_csv(args.input)
    df_bias = df[df['datatype'] == 'bias']
    df_bias_unc = df[df['datatype'] == 'bias_unc']
    df_coverage = df[df['datatype'] == 'coverage']
    df_obs = df[labels_obs]
    df_obs = df_obs.drop_duplicates(ignore_index = True)

    for _,row in df_obs.iterrows():
      bs_rows = get_bs_entries(df_bias,row)
      #print(row,row['obs1'],row['obs2'])
      bin_labels = [name for name in list(bs_rows.columns.values) if (row['obs1'] in name and row['obs2'] in name and ~np.isnan(bs_rows.iloc[0][name]))]
      fig, ax = plt.subplots(figsize=(10 + len(bin_labels)/3, 10))
      sns.boxplot(data=bs_rows[bin_labels],ax=ax,showfliers=False,whis=ONE_SIGMA_PERCS,color="cornflowerblue",labels=['quantiles'])
      bias_median = bs_rows[bin_labels].median()
      bs_rows_biasunc = get_bs_entries(df_bias_unc,row)
      bias_unc = bs_rows_biasunc[bin_labels].mean()
      ax.errorbar(np.arange(len(bin_labels))+0.2,bias_median,yerr=bias_unc,color='b',marker='o',label="average unc.",ls='none')
      ax.set_ylabel("unfold/truth")
      ax.axhline(y=1., ls="--", c="black")
      bottom,up = ax.get_ylim()
      lim = max(1.-bottom,up-1.0)
      ax.set_ylim((1-lim, (1+1.15*lim)))
      ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
      title = "Toy models unfolding / truth"
      ax.text(0.01, 0.85, title, fontsize="small", transform=ax.transAxes)
      ax.minorticks_off()
      ax.legend(facecolor='white', framealpha=0.95)
      hep.cms.label(ax=ax, data=False, label="Preliminary")
      fig.tight_layout(rect=[0, 0, 1, 0.95])
      fig.savefig(os.path.join(args.output, f"{row['obs1']}_{row['obs2']}_iter{row['iter']}_bias.pdf"))
      plt.close()

      bs_rows = get_bs_entries(df_coverage,row)
      print("len(bs_rows)",len(bs_rows))
      fig, ax = plt.subplots(figsize=(10 + len(bin_labels)/3, 10))
      sns.barplot(data=bs_rows[bin_labels],ax=ax,ci=68.2)
      x_coords = [p.get_x() + 0.5 * p.get_width() for p in ax.patches]
      y_coords = [p.get_height() for p in ax.patches]
      x_coords_adderror = []
      y_coords_adderror = []
      y_err_adderror = []
      for x,y in zip(x_coords,y_coords):
        if y==1:
          x_coords_adderror.append(x)
          low_limit = pow((1-0.682)/2,1./len(bs_rows))
          print("low_limit: ",low_limit)
          y_coords_adderror.append((1.+low_limit)/2)
          print("y_coords_adderror ",(1+low_limit)/2)
          y_err_adderror.append((1-low_limit)/2)
          print("y_err_adderror ",(1-low_limit)/2)
      ax.errorbar(x=x_coords_adderror, y=y_coords_adderror, yerr=y_err_adderror, fmt="none", c="k")
      ax.axhline(y=0.682, ls="--", c="black")
      ax.set_ylabel("coverage")
      ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
      ax.text(0.01, 0.85, title, fontsize="small", transform=ax.transAxes)
      ax.minorticks_off()
      hep.cms.label(ax=ax, data=False, label="Preliminary")
      fig.tight_layout(rect=[0, 0, 1, 0.95])
      fig.savefig(os.path.join(args.output, f"{row['obs1']}_{row['obs2']}_iter{row['iter']}_coverage.pdf"))
      plt.close()
    labels_obs = ['obs1','obs2']
    df_obs = df[labels_obs]
    df_obs = df_obs.drop_duplicates(ignore_index = True)
    for _,row in df_obs.iterrows():
      bs_rows = get_bs_entries(df_bias,row)
      bin_labels = [name for name in list(bs_rows.columns.values) if (row['obs1'] in name and row['obs2'] in name and ~np.isnan(bs_rows.iloc[0][name]))]
      bs_rows_biasunc = get_bs_entries(df_bias_unc,row)
      iters = bs_rows['iter']
      iters = iters.drop_duplicates()
      fig, ax = plt.subplots(len(bin_labels),1,figsize=(5 + len(iters)/6, 10+len(bin_labels)))
      for ibin, bin_label in enumerate(bin_labels):
        sns.boxplot(data=bs_rows[[bin_label,'iter']], x = 'iter',y = bin_label,ax=ax[ibin],showfliers=False,whis=ONE_SIGMA_PERCS,color="cornflowerblue",labels=['quantiles']).set(xlabel=None,ylabel=None)
        bias_median = [bs_rows[bin_label][bs_rows['iter']==iter].median() for iter in iters] 
        bias_unc = [bs_rows_biasunc[bin_label][bs_rows_biasunc['iter']==iter].mean() for iter in iters]
        ax[ibin].errorbar(np.arange(len(iters))+0.2,bias_median,yerr=bias_unc,color='b',marker='o',label="average unc.",ls='none')
        ax[ibin].axhline(y=1., ls="--", c="black")
        bottom,up = ax[ibin].get_ylim()
        #lim = max(1.-bottom,up-1.0)
        ax[ibin].set_ylim((bottom, up+0.08))
        title = bin_label
        ax[ibin].text(0.01, 0.80, title, fontsize="xx-small", transform=ax[ibin].transAxes)
        ax[ibin].minorticks_off()
        if ibin < len(bin_labels)-1:
          ax[ibin].set_xticklabels([])
        ax[ibin].tick_params(axis='both', which='major', labelsize=10)
      ax[-1].legend(facecolor='white', framealpha=0.95,fontsize=15,loc='upper center', bbox_to_anchor=(0.5, -0.3))
      hep.cms.label(ax=ax[0], data=False, label="Preliminary",fontsize=15)
      ax[-1].set_xlabel("iteration",fontsize=15)
      fig.text(0.05,0.5,"unfold/truth",ha='center', va='center', rotation='vertical')
      fig.tight_layout()
      fig.savefig(os.path.join(args.output, f"{row['obs1']}_{row['obs2']}_bias_evolution.pdf"))
      plt.close()

      bs_rows = get_bs_entries(df_coverage,row)
      fig, ax = plt.subplots(len(bin_labels),1,figsize=(5 + len(iters)/6, 10+len(bin_labels)))
      for ibin, bin_label in enumerate(bin_labels):
        sns.pointplot(data=bs_rows[[bin_label,'iter']],x = 'iter',y = bin_label,ax=ax[ibin],ci=68.2,color='b').set(xlabel=None,ylabel=None)
        x_coords = np.ma.getdata(ax[ibin].collections[0].get_offsets()[:, 0])
        y_coords = np.ma.getdata(ax[ibin].collections[0].get_offsets()[:, 1])
        x_coords_adderror = []
        y_coords_adderror = []
        y_err_adderror = []
        print(x_coords,y_coords)
        for x,y in zip(x_coords,y_coords):
          if y==1:
            x_coords_adderror.append(x)
            low_limit = pow((1-0.682)/2,1./(len(bs_rows)/len(iters)))
            print("low_limit: ",low_limit)
            y_coords_adderror.append((1+low_limit)/2)
            print("y_coords_adderror ",(1+low_limit)/2)
            y_err_adderror.append((1-low_limit)/2)
            print("y_err_adderror ",(1-low_limit)/2)
        ax[ibin].errorbar(x=x_coords_adderror, y=y_coords_adderror, yerr=y_err_adderror, fmt="none", c="b")
        ax[ibin].axhline(y=0.682, ls="--", c="black",label="68.2% coverage")
        bottom,up = ax[ibin].get_ylim()
        ax[ibin].set_ylim((bottom, up+0.1))
        title = bin_label
        ax[ibin].text(0.01, 0.80, title, fontsize="xx-small", transform=ax[ibin].transAxes)
        ax[ibin].minorticks_off()
        if ibin < len(bin_labels)-1:
          ax[ibin].set_xticklabels([])
        ax[ibin].tick_params(axis='both', which='major', labelsize=10)
      hep.cms.label(ax=ax[0], data=False, label="Preliminary",fontsize=15)
      ax[-1].set_xlabel("iteration",fontsize=15)
      ax[-1].legend(facecolor='white', framealpha=0.95,fontsize=15,loc='upper center', bbox_to_anchor=(0.5, -0.3))
      fig.text(0.05,0.5,"coverage",ha='center', va='center', rotation='vertical')
      fig.tight_layout()
      fig.savefig(os.path.join(args.output, f"{row['obs1']}_{row['obs2']}_coverage_evolution.pdf"))
      plt.close()
