import numpy as np
from argparse import ArgumentParser
import os
import matplotlib as mpl
mpl.use('Agg')
import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS")


def lumi_text(lumi):
    if np.log10(lumi) >= 0 and np.log10(lumi) < 3:
        lumi_print = lumi
        lumi_unit = '$\mathrm{fb}^{-1}$'
    elif np.log10(lumi) < 6 and np.log10(lumi) >= 3:
        lumi_print = lumi * 0.001
        lumi_unit = '$\mathrm{ab}^{-1}$'
    elif np.log10(lumi) < 0 and np.log10(lumi) >= -3:
        lumi_print = lumi * 1000
        lumi_unit = '$\mathrm{pb}^{-1}$'
    elif np.log10(lumi) < -3 and np.log10(lumi) >= -6:
        lumi_print = lumi * 1000000.0
        lumi_unit = '$\mathrm{nb}^{-1}$'
    elif np.log10(lumi) <= -6 and np.log10(lumi) >= -9:
        lumi_print = lumi * 1000000000.0
        lumi_unit = '$\mathrm{\mu b}^{-1}$'
    elif np.log10(lumi) <= -9 and np.log10(lumi) >= -12:
        lumi_print = lumi * 1000000000000.0
        lumi_unit = '$\mathrm{mb}^{-1}$'
    elif np.log10(lumi) <= -12 and np.log10(lumi) >= -15:
        lumi_print = lumi * 1000000000000000.0
        lumi_unit = '$\mathrm{b}^{-1}$'
    return ('{0:.1f}').format(lumi_print) + ' ' + lumi_unit + ' 2018 (13 TeV)'

def latex_root_to_mpl(text):
  text=text.replace('#','\\')
  text=text.replace(' ','\ ')
  text = "$\mathrm{" + text + "}$"
  return text
#temporary class from root histograms to arrays, will be replaced by boost hisograms
class HistArray:

  def __init__(self, histlist = None):
    self._nested_value = [np.array([hist.GetBinContent(ibin+1) for ibin in range(hist.GetNbinsX())]) for hist in histlist.root_hists] if histlist is not None else None
    self._nested_error = [np.array([hist.GetBinError(ibin+1) for ibin in range(hist.GetNbinsX())]) for hist in histlist.root_hists] if histlist is not None else None
    if  histlist is not None:
      self._nested_bins = []
      for hist in histlist.root_hists:
        bins = []
        for ibin in range(hist.GetNbinsX()):
          bins.append(hist.GetXaxis().GetBinLowEdge(ibin+1))
        bins.append(hist.GetXaxis().GetBinUpEdge(hist.GetNbinsX()))
        self._nested_bins.append(np.array(bins))
    else:
      self._nested_bins = None

  def value_array(self,i):
    return self._nested_value[i]

  def error_array(self,i):
    return self._nested_error[i]

  def bins_array(self,i):
    return self._nested_bins[i]

  def length(self):
    return len(self._nested_value)

  def __truediv__(self,other):
    result = HistArray()
    result._nested_bins = self._nested_bins
    result._nested_value = [self._nested_value[i] / other._nested_value[i] for i in range(self.length())]
    result._nested_error = [self._nested_error[i] / other._nested_value[i] for i in range(self.length())]
    return result

def MapMarkerStyle(style):
  if style == 'marker':
    return 'o'
  elif style == 'cross':
    return '.'
  elif style == 'triangle':
    return '^'
  elif style == 'triangle_down':
    return 'v'
  else:
    return None

def draw_array(value,error,bins,style,ax,color,legend):

  marker = MapMarkerStyle(style)
  if marker is not None:
    center = [(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]
    xerrs = [(bins[i+1] - bins[i])/2 for i in range(len(bins)-1)]
    ax.errorbar(
      center,
      list(value),
      xerr = xerrs,
      yerr = list(error),
      color = color,
      marker = marker,
      label = legend,
      ls = 'none'
    )
  elif style == 'fillederror':
    error_up = np.append(value,value[-1]) + np.append(error,error[-1])
    error_down = np.append(value,value[-1]) - np.append(error,error[-1])
    ax.fill_between(
      x = bins,
      y1 = error_up,
      y2 = error_down,
      color = color,
      step = "post",
      linewidth = 0,
      alpha = 0.8,
      label = legend
    )
    ax.fill_between(
      x = bins,
      y1 = np.append(value,value[-1]),
      y2 = np.append(value,value[-1]),
      color = color,
      step = "post"
    )
  elif style == 'filled':
    error_up = np.append(value,value[-1]) + np.append(error,error[-1])
    error_down = np.append(value,value[-1]) - np.append(error,error[-1])
    ax.fill_between(
      x = bins,
      y1 = error_up,
      y2 = error_down,
      color = color,
      step = "post",
      linewidth = 0,
      alpha = 0.8,
      label = legend
    )
    hep.histplot(
      value-error,
      bins=bins,
      histtype = 'fill',
      color = color,
      alpha=0.5,
      edgecolor=None,
      linewidth = 0,
      ax = ax
    )
  elif style == 'hist':
    hep.histplot(
      value,
      bins=bins,
      histtype = 'step',
      color = color,
      label = legend,
      ax = ax
    )
  elif style == 'line':
    hep.histplot(
      value,
      bins=bins,
      yerr = error,
      histtype = 'step',
      color = color,
      label = legend,
      ax = ax
    )
  else:
    print(f"style {style} not recognized")


def plot_hists(hist_ref, list_hist_compare, legend_ref, list_legend_compare, title, is_logY, do_ratio, output_path, hist_ref_stat=0, text_list=[], style_ref='marker', color_ref="black", list_style_compare=[], list_color_compare=[], labelY='Normalized Events/Bin Width', label_ratio='Data/MC',range_ratio=0.3):

  hist_ref_arrays = HistArray(hist_ref)
  hist_ref_stat_arrays = HistArray(hist_ref_stat) if hist_ref_stat != 0 else None
  list_hist_compare_arrays = [HistArray(hist_compare) for hist_compare in list_hist_compare]
  if do_ratio:
    f, axs = plt.subplots(2,hist_ref_arrays.length(),sharex=True,sharey='row',gridspec_kw={"height_ratios": (2,1) })
  else:
    f, axs = plt.subplots(1,hist_ref_arrays.length(),sharex=True,sharey='row')
  axs = axs.flatten()
  hep.cms.label('Preliminary', data=False, rlabel="", loc=2, ax = axs[0])
  hep.cms.lumitext(lumi_text(1.4817568788812e-08), ax = axs[hist_ref_arrays.length()-1])

  order_hist_array = list_hist_compare_arrays.copy()
  order_style = list_style_compare.copy()
  order_color = list_color_compare.copy()
  order_legend = list_legend_compare.copy()
  if style_ref == "filled":
     order_hist_array.insert(0,hist_ref_arrays)
     order_style.insert(0,style_ref)
     order_color.insert(0,color_ref)
     order_legend.insert(0,legend_ref)
  else:
     order_hist_array.append(hist_ref_arrays)
     order_style.append(style_ref)
     order_color.append(color_ref)
     order_legend.append(legend_ref)

  for ihist in range(hist_ref_arrays.length()):
    for (hist_array,style,color,legend) in zip(order_hist_array,order_style,order_color,order_legend):
      draw_array(hist_array.value_array(ihist), hist_array.error_array(ihist),hist_array.bins_array(ihist),style, axs[ihist],color,legend) 
      axs[ihist].text(.5,.65,latex_root_to_mpl(text_list[ihist]),horizontalalignment='center',transform=axs[ihist].transAxes)
      #if do_ratio:
      #  axs[ihist].tick_params(axis='x',labelbottom=False)
      #if ihist>0:
      #  axs[ihist].tick_params(axis='y',labelleft=False)
  axs[0].set_ylabel(labelY)
  axs[ihist].legend()
  bottom, top = axs[ihist].get_ylim()
  axs[ihist].set_ylim(top = top*1.4, bottom = 0.0)
  if do_ratio:
    order_hist_array_ratio = [hist_compare_arrays / hist_ref_arrays for hist_compare_arrays in list_hist_compare_arrays]
    order_style_ratio = ["fillederror" if style_compare == "filled" else style_compare for style_compare in list_style_compare ]
    order_color_ratio = list_color_compare.copy()
    order_hist_array_ratio.insert(0,hist_ref_arrays / hist_ref_arrays)
    order_style_ratio.insert(0,"fillederror")
    order_color_ratio.insert(0,"darkgrey")
    if hist_ref_stat_arrays is not None:
      order_hist_array_ratio.insert(1,hist_ref_stat_arrays / hist_ref_arrays)
      order_style_ratio.insert(1,"fillederror")
      order_color_ratio.insert(1,"dimgray")
    for iratio in range(hist_ref_arrays.length()):
      for (hist_array,style,color) in zip(order_hist_array_ratio,order_style_ratio,order_color_ratio):
        draw_array(hist_array.value_array(iratio), hist_array.error_array(iratio),hist_array.bins_array(iratio),style, axs[iratio+hist_ref_arrays.length()],color,None)
        #if iratio>hist_ref_arrays.length():
        #  axs[iratio].tick_params(axis='y',labelleft=False)
    axs[hist_ref_arrays.length()].set_ylabel(label_ratio)
  axs[len(axs)-1].set_xlabel(latex_root_to_mpl(title))
  axs[len(axs)-1].set_xlim(left = hist_ref_arrays.bins_array(0)[0], right = hist_ref_arrays.bins_array(0)[-1])
  ratio_bottom, ratio_top = axs[len(axs)-1].get_ylim()
  ratio_range_sym = max(1.0-ratio_bottom,ratio_top-1.0)
  ratio_bottom_sym, ratio_top_sym = 1.0-ratio_range_sym, 1.0+ratio_range_sym
  axs[len(axs)-1].set_ylim(bottom = max(0,ratio_bottom_sym), top = ratio_top_sym)
  plt.subplots_adjust(wspace=0, hspace=0)
  plt.savefig(output_path + '.pdf')
  plt.savefig(output_path + '.png')

