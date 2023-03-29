from scipy.stats import chi2 as CHI2
import numpy as np
import itertools


def check_shape(nested_list1,nested_list2):
  assert len(nested_list1) == len(nested_list2)
  for (list1,list2) in zip(nested_list1,nested_list2):
    assert len(list1) == len(list2)


def GOF_binned_from_np(values_1,values_2,errors_1=None,errors_2=None):

  if errors_1 is None and errors_2 is None:
    raise ValueError("errors_1 and errors_2 cannot be both None")

  check_shape(values_1,values_2)
  if errors_1 is not None:
    check_shape(values_1,errors_1)
  if errors_2 is not None:
    check_shape(values_2,errors_2)

  values_1 = np.array(list(itertools.chain(*values_1)))
  values_2 = np.array(list(itertools.chain(*values_2)))

  if errors_1 is None:
    errors_1=np.zeros(shape=np.shape(values_1))
  else:
    errors_1 = np.array(list(itertools.chain(*errors_1)))

  if errors_2 is None:
    errors_2=np.zeros(shape=np.shape(values_2))
  else:
    errors_2 = np.array(list(itertools.chain(*errors_2)))

  diff = values_1 - values_2
  error_square = np.square(errors_1) + np.square(errors_2)
  return GOF(diff,error_square)

def GOF_binned_from_root(hist_list_1, hist_list_2, use_error='both'):
  
  if use_error not in ['both','first','second']:
    raise ValueError("use_error can only be 'both','first' or 'second'")
  values_1 = []
  values_2 = []
  error_square = []
  for (hist1,hist2) in zip(hist_list_1,hist_list_2):
    assert hist1.GetNbinsX() == hist2.GetNbinsX()
    values_1 += [hist1.GetBinContent(ibin+1) for ibin in range(hist1.GetNbinsX())]
    values_2 += [hist2.GetBinContent(ibin+1) for ibin in range(hist2.GetNbinsX())]
    if use_error == 'first':
      error_square += [hist1.GetBinError(ibin+1)**2 for ibin in range(hist1.GetNbinsX())]
    elif use_error == 'second':
      error_square += [hist2.GetBinError(ibin+1)**2 for ibin in range(hist2.GetNbinsX())]
    else:
      error_square += [hist1.GetBinError(ibin+1)**2+hist2.GetBinError(ibin+1)**2 for ibin in range(hist1.GetNbinsX())]

  values_1 = np.array(values_1)
  values_2 = np.array(values_2)
  diff =  values_1-values_2
  error_square = np.array(error_square)
  return GOF(diff,error_square)


def GOF(diff,error_square):

  wchi2 = np.sum(diff*diff/error_square)
  ndof = len(diff)-1

  p = CHI2.sf(wchi2,ndof)
  return wchi2,ndof,wchi2/ndof,p









