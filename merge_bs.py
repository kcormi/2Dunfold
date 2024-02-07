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

def write_mean_std(entry_merge,entries_origin,entry_nominal=None):
  bin_values = entries_origin['bin_values']
  bin_values = df_column_to_ak(bin_values)
  dim1_underflow = ak.Array(list(entries_origin['dim1_underflow'].values))
  dim1_overflow = ak.Array(list(entries_origin['dim1_overflow'].values))
  if entry_nominal is not None:
    bin_values_nominal = entry_nominal['bin_values']
    bin_values_nominal = df_column_to_ak(bin_values_nominal)
    dim1_underflow_nominal = np.array(entry_nominal['dim1_underflow'].values)
    dim1_overflow_nominal = np.array(entry_nominal['dim1_overflow'].values)
    sum_nominal = np.sum(bin_values_nominal[0])+dim1_underflow_nominal[0]+dim1_overflow_nominal[0]
    sum_bs = ak.sum(ak.flatten(bin_values,axis=-1),axis=-1)+dim1_underflow+dim1_overflow
    bin_values = bin_values/sum_bs*sum_nominal
    dim1_underflow =  dim1_underflow/sum_bs*sum_nominal
    dim1_overflow =  dim1_overflow/sum_bs*sum_nominal
  entry_merge['bin_values'] = ak.to_list(ak.mean(bin_values,axis = 0))
  print(entry_merge['datatype'],entry_merge['obs1'],entry_merge['obs2'],"iter",entry_merge['iter'])
  print("bin values",bin_values)
  print("ak.mean(bin_values*bin_values,axis = 0)",*ak.mean(bin_values*bin_values,axis = 0))
  print("np.square(ak.mean(bin_values,axis = 0))",*np.square(ak.mean(bin_values,axis = 0)))
  entry_merge['bin_errors'] = ak.to_list(np.sqrt(np.absolute((ak.mean(bin_values*bin_values,axis = 0) - np.square(ak.mean(bin_values,axis = 0))))*len(bin_values)/(len(bin_values)-1)))
  print("entry_merge['bin_errors']",*entry_merge['bin_errors'])
  entry_merge['dim1_underflow'] = ak.to_list(ak.mean(dim1_underflow,axis = 0))
  entry_merge['dim1_overflow'] = ak.to_list(ak.mean(dim1_overflow,axis = 0))

def write_sum_std(entry_merge,entries_origin,entry_nominal):
  bin_values = entries_origin['bin_values']
  bin_values = df_column_to_ak(bin_values)
  dim1_underflow = ak.Array(list(entries_origin['dim1_underflow'].values))
  dim1_overflow = ak.Array(list(entries_origin['dim1_overflow'].values))
  bin_values_nominal = entry_nominal['bin_values']
  bin_values_nominal = df_column_to_ak(bin_values_nominal)
  bin_errors_nominal = entry_nominal['bin_errors']
  bin_errors_nominal = df_column_to_ak(bin_errors_nominal)
  dim1_underflow_nominal = np.array(entry_nominal['dim1_underflow'].values)
  dim1_overflow_nominal = np.array(entry_nominal['dim1_overflow'].values)
  entry_merge['bin_errors'] = ak.to_list(np.sqrt(np.square(bin_errors_nominal)[0]+ak.sum(np.square(bin_values-bin_values_nominal),axis = 0)))

def write_sum_error(entry_merge,entry_1,entry_2):
  entry_merge['bin_errors'] = ak.to_list(np.sqrt(np.square(ak.Array(entry_1['bin_errors']))+np.square(ak.Array(entry_2['bin_errors']))))



def write_cov_corr(entry_merge_cov,entry_merge_corr,entries_origin,entry_nominal_1=None,entry_nominal_2=None):
  bin_values_x = entries_origin['bin_values_x']
  bin_values_y = entries_origin['bin_values_y']
  bin_values_x = df_column_to_ak(bin_values_x)
  bin_values_y = df_column_to_ak(bin_values_y)
  if (entry_nominal_1 is not None) and (entry_nominal_2 is not None):
    bin_values_x_nominal = entry_nominal_1['bin_values']
    bin_values_y_nominal = entry_nominal_2['bin_values']
    bin_values_x_nominal = df_column_to_ak(bin_values_x_nominal)
    bin_values_y_nominal = df_column_to_ak(bin_values_y_nominal)
    sum_x_nominal = np.sum(bin_values_x_nominal[0])
    sum_y_nominal = np.sum(bin_values_y_nominal[0])
    sum_x = ak.sum(ak.flatten(bin_values_x,axis=-1),axis=-1)
    sum_y = ak.sum(ak.flatten(bin_values_y,axis=-1),axis=-1)
    bin_values_x = bin_values_x/sum_x*sum_x_nominal
    bin_values_y = bin_values_y/sum_y*sum_y_nominal
  bin_values_x =   bin_values_x[:,:,1:-1] #remove the underflow and overflow part in dim2
  bin_values_y =   bin_values_y[:,:,1:-1]
  bin_values_x = ak.flatten(bin_values_x,axis=-1)
  bin_values_y = ak.flatten(bin_values_y,axis=-1)
  if (entry_nominal_1 is not None) and (entry_nominal_2 is not None):
    sum_x = ak.sum(bin_values_x,axis=-1,keepdims=True)
    sum_y = ak.sum(bin_values_y,axis=-1,keepdims=True)
    bin_values_x = ak.concatenate([sum_x,bin_values_x],axis=1)
    bin_values_y = ak.concatenate([sum_y,bin_values_y],axis=1)
  cov = np.cov(bin_values_x,bin_values_y,rowvar=False)[:np.shape(bin_values_x)[-1],np.shape(bin_values_x)[-1]:]
  entry_merge_cov['bin_values'] = cov.to_list()
  entry_merge_cov['bin_errors'] = None
  corr = np.nan_to_num(np.corrcoef(bin_values_x,bin_values_y,rowvar=False))[:np.shape(bin_values_x)[-1],np.shape(bin_values_x)[-1]:]
  entry_merge_corr['bin_values'] = corr.to_list()
  entry_merge_corr['bin_errors'] = None

if __name__=="__main__":

    parser = ArgumentParser()
    parser.add_argument('-i','--input', default="config/merge_unfolddata.txt", help="The text file including the names of the csv files for merging")
    parser.add_argument('-o','--output', default="merge.csv", help="The name of the csv file for output")
    args = parser.parse_args()

    with open(args.input) as file:
      lines = [line.rstrip() for line in file]

    labels = ['obs1','obs2','histtype','iter','dataset','method','from_step1','do_eff_acc','dim1_isgen','dim2_isgen','bin_edges_dim1','bin_edges_dim2']
    labels_unfold = ['obs1','obs2','iter','method','from_step1','do_eff_acc','bin_edges_dim1','bin_edges_dim2']
    MCname = "A3P"


    df_merge = None
    df_sysbs_merge = None
    df_sys_unfold_merge = None
    df_statbs_merge = None

    for name in lines:
      df = pd.read_csv(name)
      df_non_bs = df[~df['datatype'].str.contains('bs')]
      df_sysbs = df[df['datatype'].str.contains('sysbs')]
      df_statbs = df[df['datatype'].str.contains('statbs')]
      df_sys_unfold = df[(df['datatype']=='unfold')&(df['histtype']=='inclusive')&df['dim1_isgen']&df['dim2_isgen']]
      df_merge = append_df(df_merge,df_non_bs)
      df_sysbs_merge = append_df(df_sysbs_merge,df_sysbs)
      df_sys_unfold_merge = append_df(df_sys_unfold_merge,df_sys_unfold)
      df_statbs_merge = append_df(df_statbs_merge,df_statbs)


    df_sys_unfold_labels = df_sys_unfold_merge[labels_unfold]
    df_sys_unfold_labels = df_sys_unfold_labels.drop_duplicates(ignore_index = True)
    for _, row in df_sys_unfold_labels.iterrows():
      print("row: ",row)
      unfold_rows = get_bs_entries(df_sys_unfold,row)
      nominal_row = unfold_rows[unfold_rows['dataset']==MCname]
      sys_rows = unfold_rows[unfold_rows['dataset']!=MCname]
      print("nominal_row: ",nominal_row)
      print("sys_rows: ",sys_rows)
      nominal_row_merge = copy.deepcopy(nominal_row.iloc[0])
      nominal_row_merge['datatype'] = 'unfold_fullunc'
      if not sys_rows.empty:
        write_sum_std(nominal_row_merge,sys_rows,nominal_row)
        df_merge = df_merge.append(copy.deepcopy(nominal_row_merge),ignore_index = True)

    for datatype,df_bs_merge in zip(['unfold_sys','unfold_stat'],[df_sysbs_merge,df_statbs_merge]):
      df_labels = df_bs_merge[labels]
      df_labels = df_labels.drop_duplicates(ignore_index = True)
      df_labels_cov = df_labels[(df_labels['histtype'] == 'inclusive')&(df_labels['dim1_isgen'])&(df_labels['dim2_isgen'])]
      for _, row in df_labels.iterrows():
        bs_rows = get_bs_entries(df_bs_merge,row)
        nominal_rows = get_bs_entries(df_merge[(df_merge['datatype']=='unfold')],row)
        bs_row_merge = copy.deepcopy(bs_rows.iloc[0])
        bs_row_merge['datatype'] = datatype
        write_mean_std(bs_row_merge,bs_rows)
        df_merge = df_merge.append(copy.deepcopy(bs_row_merge),ignore_index = True)
        bs_row_merge['datatype'] = datatype+"_fixXS"
        write_mean_std(bs_row_merge,bs_rows,entry_nominal=nominal_rows)
        df_merge = df_merge.append(copy.deepcopy(bs_row_merge),ignore_index = True)
      for rows in combinations_with_replacement(df_labels_cov.iterrows(),2):
        _,row1 = rows[0]
        _,row2 = rows[1]
        if row1['iter'] != row2['iter']: continue
        bs_rows_1 = get_bs_entries(df_bs_merge,row1)
        bs_rows_2 = get_bs_entries(df_bs_merge,row2)
        nominal_row_1 = get_bs_entries(df_merge[(df_merge['datatype']=='unfold')],row1)
        nominal_row_2 = get_bs_entries(df_merge[(df_merge['datatype']=='unfold')],row2)
        cov_row  = copy.deepcopy(bs_rows_1.iloc[0])
        for column_obs in ['obs1','obs2','bin_edges_dim1','bin_edges_dim2']:
          cov_row[column_obs] = (row1[column_obs],row2[column_obs])
        cov_row['datatype'] = datatype
        corr_row = copy.deepcopy(cov_row)
        cov_row['histtype'] = 'covariance'
        corr_row['histtype'] = 'correlation'
        bs_rows_combine = bs_rows_1.merge(bs_rows_2,on='datatype')
        write_cov_corr(cov_row,corr_row,bs_rows_combine)
        df_merge = df_merge.append(copy.deepcopy(cov_row),ignore_index = True)
        df_merge = df_merge.append(copy.deepcopy(corr_row),ignore_index = True)
        cov_row['histtype'] = 'covariance_fixXS'
        corr_row['histtype'] = 'correlation_fixXS'
        write_cov_corr(cov_row,corr_row,bs_rows_combine,entry_nominal_1=nominal_row_1,entry_nominal_2=nominal_row_2)
        df_merge = df_merge.append(copy.deepcopy(cov_row),ignore_index = True)
        df_merge = df_merge.append(copy.deepcopy(corr_row),ignore_index = True)

    if df_merge[df_merge['datatype']=='unfold_stat_fixXS'].size>0 and df_merge[df_merge['datatype']=='unfold_sys_fixXS'].size>0:
      for _, row in df_sys_unfold_labels.iterrows():
        unfold_rows = get_bs_entries(df_sys_unfold,row)
        nominal_row = unfold_rows[unfold_rows['dataset']==MCname]
        bs_row_mergestatsys = copy.deepcopy(nominal_row.iloc[0])
        bs_row_mergestatsys['datatype'] = 'unfold_fullbsunc_fixXS'
        statbsrow_merge = get_bs_entries(df_merge[df_merge['datatype']=='unfold_stat_fixXS'],row)
        sysbsrow_merge = get_bs_entries(df_merge[df_merge['datatype']=='unfold_sys_fixXS'],row) 
        write_sum_error(bs_row_mergestatsys,statbsrow_merge.iloc[0],sysbsrow_merge.iloc[0])
        df_merge = df_merge.append(copy.deepcopy(bs_row_mergestatsys))
    df_merge.to_csv(args.output,mode='w',index=False)

