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

def append_df(df_origin,df_append):
  if df_origin is None:
    df_origin = copy.deepcopy(df_append)
  else:
    df_origin = df_origin.append(df_append,ignore_index = True)
  df_origin.drop_duplicates(ignore_index = True, inplace = True)
  return df_origin

def get_bs_entries(df,label_row):
  for ikey, key in enumerate(list(label_row.keys())):
    if ikey == 0:
      sel = (df[key]==label_row[key])
    else:
      sel = sel & (df[key]==label_row[key])
  return df[sel]

def df_column_to_ak(column):
  array = column.replace('inf','2e308')
  array = array.apply(ast.literal_eval)
  array = array.values
  array = list(array)
  array = ak.Array(array)
  return array

def write_bias_coverage(entry,entry_bias,entry_bias_unc,entry_coverage):
  bin_edges_dim1 = entry['bin_edges_dim1']
  bin_edges_dim1 = ast.literal_eval(bin_edges_dim1)
  bin_edges_dim2 = entry['bin_edges_dim2']
  bin_edges_dim2 = ast.literal_eval(bin_edges_dim2)
  obs1 = entry['obs1']
  obs2 = entry['obs2']
  flat_bin_list = []
  bias_list = []
  bias_unc_list = []
  coverage_list = []
  for ibin_dim1,lowbin_dim1 in enumerate(bin_edges_dim1[:-1]):
    for ibin_dim2,lowbin_dim2 in enumerate(bin_edges_dim2[ibin_dim1][:-1]):
      flat_bin_list.append(f'{obs1} {lowbin_dim1}-{bin_edges_dim1[ibin_dim1+1]} {obs2} {lowbin_dim2}-{bin_edges_dim2[ibin_dim1][ibin_dim2+1]}')
      bias_list.append(entry['bin_values'][ibin_dim1][ibin_dim2])
      bias_unc_list.append(entry['bin_errors'][ibin_dim1][ibin_dim2])
      coverage_list.append(int(math.fabs(entry['bin_values'][ibin_dim1][ibin_dim2]-1)<entry['bin_errors'][ibin_dim1][ibin_dim2]))
  for flat_bin, bias, bias_unc, coverage in zip(flat_bin_list,bias_list,bias_unc_list,coverage_list):
    entry_bias[flat_bin] = bias
    entry_bias_unc[flat_bin] = bias_unc
    entry_coverage[flat_bin] = coverage



def write_ratio(entry_ratio,row_truth,row_unfold,row_unc):
  bin_values_truth = row_truth['bin_values']
  bin_values_truth = df_column_to_ak(bin_values_truth)
  bin_values_unfold = row_unfold['bin_values']
  bin_values_unfold = df_column_to_ak(bin_values_unfold)
  bin_errors_unfold = row_unc['bin_errors']
  bin_errors_unfold = df_column_to_ak(bin_errors_unfold)
  entry_ratio['bin_values'] = ak.to_list(bin_values_unfold[0][:,1:-1]/bin_values_truth[0][:,1:-1])
  entry_ratio['bin_errors'] = ak.to_list(bin_errors_unfold[0][:,1:-1]/bin_values_truth[0][:,1:-1])



if __name__=="__main__":

    parser = ArgumentParser()
    parser.add_argument('-i','--input', default="config/merge_unfoldCP1_stat_sys_bs.txt", help="The text file including the names of the csv files for merging")
    parser.add_argument('-unc','--uncertainty', default="results_finebin_v7_MCA3P_unfolddata_v3_2/merge.csv", help="The text file including the names of the csv files for merging")
    parser.add_argument('--dataunc',action="store_true",help="Use the uncertainty from data unfolding with bootstraps")
    parser.add_argument('-o','--output', default="results_finebin_v7_MCEPOS_unfoldCP1_1d_eff_acc_ensemble4_coverage_v2/merge_bs.csv", help="The name of the csv file for output")
    args = parser.parse_args()

    if args.dataunc:
      if os.path.isfile(args.uncertainty):
        try:
          df_unc = pd.read_csv(args.uncertainty)
          df_unc = df_unc[df_unc['datatype']=='unfold_fullbsunc_fixXS']
        except ValueError:
          print("file "+args.uncertainty+" does not exist")

    with open(args.input) as file:
      lines = [line.rstrip() for line in file]

    df_merge = None
    df_labels_merge = None


    labels = ['obs1','obs2','iter','dataset','method','from_step1','do_eff_acc']
    labels_pseudodata = ['obs1','obs2','method','from_step1','do_eff_acc']
    diclist_bias_coverage = []

    for name in lines:
      df = pd.read_csv(name)
      df_pseudodata = df[(df["datatype"]=="pseudodata")&(df['histtype']=='inclusive')&df['dim1_isgen']&df['dim2_isgen']]
      df_unfold = df[df['datatype']=='unfold_fullunc']
      df_labels = df_unfold[labels]
      df_labels = df_labels.drop_duplicates(ignore_index = True)
      df_labels_merge = append_df(df_labels_merge,df_labels)
      for _, row in df_labels.iterrows():
        row_pseudodata = get_bs_entries(df_pseudodata,row[labels_pseudodata])
        row_unfold = get_bs_entries(df_unfold,row)
        if args.dataunc:
          row_unc = get_bs_entries(df_unc,row)
        else:
          row_unc = row_unfold
        entry_ratio = copy.deepcopy(row_unfold.iloc[0])
        entry_ratio["datatype"] = "unfold_to_truth"
        write_ratio(entry_ratio,row_pseudodata,row_unfold,row_unc)
        entry_bias = copy.deepcopy(entry_ratio)
        entry_bias["datatype"] = "bias"
        entry_bias_unc = copy.deepcopy(entry_ratio)
        entry_bias_unc["datatype"] = "bias_unc"
        entry_coverage = copy.deepcopy(entry_ratio)
        entry_coverage["datatype"] = "coverage"
        write_bias_coverage(entry_ratio,entry_bias,entry_bias_unc,entry_coverage)
        diclist_bias_coverage.append(entry_bias)
        diclist_bias_coverage.append(entry_bias_unc)
        diclist_bias_coverage.append(entry_coverage)
    pd.DataFrame(diclist_bias_coverage).to_csv(args.output,mode='w',index=False)

