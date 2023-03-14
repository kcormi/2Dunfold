# uncompyle6 version 3.9.0
# Python bytecode version base 2.7 (62211)
# Decompiled from: Python 3.7.1 (default, Dec 14 2018, 13:28:58) 
# [Clang 4.0.1 (tags/RELEASE_401/final)]
# Embedded file name: /work/jinw/CMSSW_10_2_15_patch2/src/cleanup_2Dunfold/unfold_utils.py
# Compiled at: 2023-03-06 20:32:07
import numpy as np, json
from argparse import ArgumentParser
import os 
import ROOT 
from enum import Enum

class CutType(Enum):
    PassReco = 1
    PassGen  = 2 
    PassReco_PassGen = 3
    PassReco_noPassGen = 4
    noPassReco_PassGen = 5
    noPassReco_noPassGen = 6

class HistType(Enum):
    RECO = 1
    GEN = 2
    MIG = 3

root_cuts = {}
root_cuts[CutType.PassReco] = '(PV_N_good==1&&PV_isgood&&(Instanton_N_Trk_highPurity_pt05>2))'
root_cuts[CutType.PassGen] = '(Instanton_N_gen_ChargedFSParticle_eta2p4pt05>2)'
root_cuts[CutType.PassReco_PassGen] = '(PV_N_good==1&&PV_isgood&&(Instanton_N_Trk_highPurity_pt05>2)&&(Instanton_N_gen_ChargedFSParticle_eta2p4pt05>2))'
root_cuts[CutType.PassReco_noPassGen] = '(PV_N_good==1&&PV_isgood&&(Instanton_N_Trk_highPurity_pt05>2)&&(Instanton_N_gen_ChargedFSParticle_eta2p4pt05<=2))'
root_cuts[CutType.noPassReco_PassGen] = '(((PV_N_good!=1)||(PV_isgood&&(Instanton_N_Trk_highPurity_pt05<=2)))&&(Instanton_N_gen_ChargedFSParticle_eta2p4pt05>2))'
root_cuts[CutType.noPassReco_noPassGen] = '(((PV_N_good!=1)||(PV_isgood&&(Instanton_N_Trk_highPurity_pt05<=2)))&&(Instanton_N_gen_ChargedFSParticle_eta2p4pt05<=2))'

np_cuts = {}
np_cuts[CutType.PassReco] = [['reco_ntrk', '!=', 'NAN']]
np_cuts[CutType.PassGen] = [['gen_nch', '>', '2']]
np_cuts[CutType.PassReco_PassGen] = [['reco_ntrk', '!=', 'NAN'], ['gen_nch', '>', '2']]
np_cuts[CutType.PassReco_noPassGen] = [['reco_ntrk', '!=', 'NAN'], ['gen_nch', '<=', '2']]
np_cuts[CutType.noPassReco_PassGen] = [['reco_ntrk', '==', 'NAN'], ['gen_nch', '>', '2']]
np_cuts[CutType.noPassReco_noPassGen] = [['reco_ntrk', '==', 'NAN'], ['gen_nch', '<=', '2']]

def filter_np_cut(np_dic, np_cut):
    filter_cut = np.ones(len(np_dic['reco_ntrk']), dtype=bool)
    if np_cut is None:
        return filter_cut
    else:
        for cut in np_cut:
            if cut[1] == '!=':
                if np.isnan(float(cut[2])):
                    filter_cut = filter_cut * (np.isnan(np_dic[cut[0]]) == False)
                else:
                    filter_cut = filter_cut * (np_dic[cut[0]] != float(cut[2]))
            elif cut[1] == '==':
                if np.isnan(float(cut[2])):
                    filter_cut = filter_cut * np.isnan(np_dic[cut[0]])
                else:
                    filter_cut = filter_cut * (np_dic[cut[0]] == float(cut[2]))
            elif cut[1] == '>':
                filter_cut = filter_cut * (np_dic[cut[0]] > float(cut[2]))
            elif cut[1] == '>=':
                filter_cut = filter_cut * (np_dic[cut[0]] >= float(cut[2]))
            elif cut[1] == '<':
                filter_cut = filter_cut * (np_dic[cut[0]] < float(cut[2]))
            elif cut[1] == '<=':
                filter_cut = filter_cut * (np_dic[cut[0]] <= float(cut[2]))
            else:
                raise ValueError('Cannot recognize the cut ', cut)

        return filter_cut


def merge_bins(obs, trees=[], root_cut='', threshold=1.0, bin_edges_dim1_1d=None, bin_edges_dim2_1d=None):
    bin_edges_dim1_1d = [] if bin_edges_dim1_1d is None else bin_edges_dim1_1d
    bin_edges_dim2_1d = [] if bin_edges_dim2_1d is None else bin_edges_dim2_1d
    hists = []
    for itree, tree in enumerate(trees):
        hist = ROOT.TH2F(f'HistMerge{itree}', f'HistMerge{itree}', len(bin_edges_dim1_1d) - 1, np.array(bin_edges_dim1_1d), len(bin_edges_dim2_1d) - 1, np.array(bin_edges_dim2_1d))
        tree.Draw(f'{obs[1]}:{obs[0]}>>HistMerge{itree}', root_cut, 'colzgoff')
        hists.append(hist)

    bin_edges_dim2_merge = []
    for ibin0 in range(len(bin_edges_dim1_1d) - 1):
        bin_edges_dim2_merge.append(list(bin_edges_dim2_1d))
        ibin1 = 0
        content = np.zeros(len(hists))
        merge_flag = True
        while ibin1 < len(bin_edges_dim2_1d) - 1:
            content = content + np.array([ hist.GetBinContent(ibin0 + 1, ibin1 + 1) for hist in hists ])
            if np.all(content > threshold):
                content = np.zeros(len(hists))
                merge_flag = False
            else:
                merge_flag = True
                if ibin1 != len(bin_edges_dim2_1d) - 2:
                    bin_edges_dim2_merge[ibin0].remove(bin_edges_dim2_1d[ibin1 + 1])
            ibin1 += 1

        if np.any(content <= threshold) and merge_flag:
            if len(bin_edges_dim2_merge[ibin0]) > 2:
                del bin_edges_dim2_merge[ibin0][-2]
            else:
                raise ValueError('Not enough statistics in first observable bin ' + str(ibin0))

    return bin_edges_dim2_merge


class Cut:

    def __init__(self, root_cut, npy_cut=None):
        self.root_cut = root_cut
        self.npy_cut = npy_cut

cuts = {}
for cut_type in CutType:
    cuts[cut_type] = Cut( root_cuts[cut_type], np_cuts[cut_type])


class HistDim:
    
    def __init__(self, is_gen=False, inner_dim=None):
        self.is_gen = is_gen
        self.bin_edges = None
        self.root_var = ''
        self.np_var  = ''
        self.underflow = 0.
        self.overflow = 0.
        self.inner_dim = inner_dim

    def from_config(self, config, isgen=False, inner_dim = None):
        inner_dim_factor = 0
        self.inner_dim = inner_dim
        if self.inner_dim is not None:
            inner_dim_factor = max(1, self.inner_dim.length() )
        if isgen:
            if self.inner_dim is None:
                self.bin_edges = config['binedgesgen'] 
            else:
                self.bin_edges = [config['binedgesgen']] * inner_dim_factor
            self.root_var = config['gen']
            self.np_var = config['gen_key']
        else:
            if self.inner_dim is None:
                self.bin_edges = config['binedgesreco'] 
            else:
                self.bin_edges = [config['binedgesreco']] * inner_dim_factor
            self.root_var = config['reco']
            self.np_var = config['reco_key']

    def length(self):
        if self.bin_edges is None or len(self.bin_edges) == 0: 
            return 0
        return len(self.bin_edges) - 1



class HistList:


    def __init__(self, name, tag=''):
        self.htype = None
        self._name = name
        self.tag = tag
        self.legend = ''
        self.color = ROOT.kBlack
        self.style = None
        self.bin_values = None
        self.dim1 = HistDim()
        self.dim2 = HistDim()
        self.cut = Cut('',None)
        self.root_hists = None
        self.root_2Dhist = None
        self.root_flat_hist = None
        self.flat_xAxisLabel = None
        self.flat_xAxisLabel_up = None
        self.flat_xAxisLabel_low = None
        self.root_hists_name = []
        self.bin_sum = []
        self.bin_norm = []
        self.total = 0.0
        self.dim1_underflow = 0.0
        self.dim1_overflow = 0.0
        self.norm = 0.0
        return

    @property
    def name(self):
        return self._name

    @property
    def nd1(self):
        return self.dim1.length()

    @property
    def nd2(self):
        return self.dim1.length()

    @property
    def root_cut(self):
        return self.cut.root_cut

    @root_cut.setter
    def root_cut(self, root_cut):
        self.cut.root_cut = root_cut

    @property
    def npy_cut(self):
        return self.cut.npy_cut

    @npy_cut.setter
    def npy_cut(self, npy_cut):
        self.cut.npy_cut = npy_cut

    def fill_root_hists_name(self):
        for i in range(self.dim1.length()):
            self.root_hists_name.append(f'{self.name}_{self.dim1.np_var}_{i+1}{self.tag}')

    @property
    def bin_edges_dim1(self):
        return self.dim1.bin_edges

    @bin_edges_dim1.setter
    def bin_edges_dim1(self, edges):
        self.dim1.bin_edges = edges

    @property
    def bin_edges_dim2(self):
        return self.dim2.bin_edges

    @bin_edges_dim2.setter
    def bin_edges_dim2(self, edges):
        self.dim2.bin_edges = edges

    @property
    def root_variable_dim1(self):
        return self.dim1.root_var

    @property
    def root_variable_dim2(self):
        return self.dim2.root_var

    @property
    def np_variable_dim1(self):
        return self.dim1.np_variable

    @property
    def np_variable_dim2(self):
        return self.dim2.np_variable


    def read_settings_from_config_dim1(self, config, isgen=False):
        self.dim1.from_config( config, isgen=isgen )

    def read_settings_from_config_dim2(self, config, isgen=False):
        self.dim2.from_config( config, isgen=isgen, inner_dim=self.dim1 )

    def fill_hist( self, event_info, from_root, weightarray=None, genWeight=''):
        if from_root:
            self.fill_hist_from_root( event_info, genWeight=genWeight)
        else:
            self.fill_hist_from_npz( event_info, weightarray= weightarray, genWeight=genWeight )

    def fill_hist_from_root(self, tree, genWeight=''):
        if len(self.dim1.bin_edges) - 1 != len(self.dim2.bin_edges):
            raise ValueError(('bin_edges_dim1 (length {0}) and bin_edges_dim2 (length {1}) not compatible').format(len(self.dim1.bin_edges), len(self.dim2.bin_edges)))
        if len(self.root_hists_name) != len(self.dim2.bin_edges):
            raise ValueError('root_hists_name does not have the same length as bin_edges_dim2. Run fill_root_hists_name() first')
        self.root_hists = [ ROOT.TH1F(hist_name, hist_name, len(self.dim2.bin_edges[ibin1]) - 1, np.array(self.dim2.bin_edges[ibin1])) for ibin1, hist_name in enumerate(self.root_hists_name) ]
        string_weight = ''
        if genWeight != '':
            string_weight = '*' + genWeight
        base_cut_str = self.root_cut + string_weight
        for ihist, edge_i in enumerate(self.dim1.bin_edges[:-1]): 
            edge_ip1 = self.dim1.bin_edges[ihist + 1]
            cut_str = f'{base_cut_str}*({self.dim1.root_var}>={edge_i})*({self.dim1.root_var}<{edge_ip1})'
            print(cut_str)
            print(tree.GetEntries(cut_str))
            tree.Draw(f'{self.dim2.root_var}>>{self.root_hists_name[ihist]}', cut_str)
            self.bin_sum.append(np.sum([ self.root_hists[ihist].GetBinContent(ibin + 1) for ibin in range(self.root_hists[ihist].GetNbinsX()) ]))
            self.bin_norm.append(np.sum([ self.root_hists[ihist].GetBinContent(ibin) for ibin in range(self.root_hists[ihist].GetNbinsX() + 2) ]))

        self.dim1_underflow = tree.GetEntries(f'{base_cut_str}*({self.dim1.root_var}<{self.dim1.bin_edges[0]})')
        self.dim1_overflow = tree.GetEntries(f'{base_cut_str}*({self.dim1.root_var}>={self.dim1.bin_edges[-1]})')
        self.total = np.sum(self.bin_sum)
        self.norm = np.sum(self.bin_norm) + self.dim1_overflow + self.dim1_underflow

    def fill_hist_from_npz(self, files, weightarray=None, genWeight=''):
        if len(self.dim1.bin_edges) - 1 != len(self.dim2.bin_edges):
            raise ValueError(('bin_edges_dim1 (length {0}) and bin_edges_dim2 (length {1}) not compatible').format(len(self.dim1.bin_edges), len(self.dim2.bin_edges)))
        if len(self.root_hists_name) != len(self.dim2.bin_edges):
            raise ValueError('root_hists_name does not have the same length as bin_edges_dim2. Run fill_root_hists_name() first')
        self.root_hists = [ ROOT.TH1F(hist_name, hist_name, len(self.dim2.bin_edges[ibin1]) - 1, np.array(self.dim2.bin_edges[ibin1])) for ibin1, hist_name in enumerate(self.root_hists_name) ]
        if isinstance(files, list):
            f_list = [ np.load(str(file_one), allow_pickle=True) for file_one in files ]
            obs_arrays = {}
            for key in list(f_list[0].keys()):
                if key != 'tracks' and key != 'charged':
                    obs_arrays[key] = np.concatenate([ f[key] for f in f_list ], axis=0)

        else:
            obs_arrays = np.load(files)
        if weightarray is None:
            weightarray = np.ones(len(obs_arrays['reco_ntrk']))
        if genWeight != '':
            if genWeight in obs_arrays.keys():
                weight = weightarray * obs_arrays[genWeight]
            else:
                print(f'{genWeight} is not a key in {files} skip it')
                weight = weightarray
        else:
            weight = weightarray
        filter_cut = filter_np_cut(obs_arrays, self.npy_cut)
        inf_weight_mask = np.isinf(weight) != True

        for ihist, edge_i in enumerate(self.dim1.bin_edges[:-1]): 
            edge_ip1 = self.dim1.bin_edges[ihist+1]
            bind1_low_cut = obs_arrays[self.dim1.np_var] >= edge_i 
            bind1_high_cut = obs_arrays[self.dim1.np_var] < edge_ip1 

            bd1_cut = bind1_low_cut & bind1_high_cut
            d1_cut = filter_cut & inf_weight_mask & bd1_cut

            for ibin in range(self.root_hists[ihist].GetNbinsX()):
                bind2_low_cut = obs_arrays[self.dim2.np_var] >= self.dim2.bin_edges[ihist][ibin]
                bind2_high_cut = obs_arrays[self.dim2.np_var] < self.dim2.bin_edges[ihist][ibin+1]
                bd2_cut = bind2_low_cut & bind2_high_cut

                full_sel = d1_cut & bd2_cut

                self.root_hists[ihist].SetBinContent(ibin + 1, np.sum(weight[full_sel]) )
                self.root_hists[ihist].SetBinError( ibin + 1, np.sqrt(np.sum(np.square(weight[full_sel]))) )

            d2_undfl_cut = obs_arrays[self.dim2.np_var] < self.dim2.bin_edges[ihist][0]
            full_sel = d1_cut & d2_undfl_cut

            self.root_hists[ihist].SetBinContent(0, np.sum(weight[full_sel]) )
            self.root_hists[ihist].SetBinContent(0, np.sqrt(np.sum(np.square(weight[full_sel]))) )

            d2_ovfl_cut = obs_arrays[self.dim2.np_var ] > self.dim2.bin_edges[ihist][-1]
            full_sel = d1_cut & d2_ovfl_cut

            self.root_hists[ihist].SetBinContent(self.root_hists[ihist].GetNbinsX() + 1, np.sum(weight[full_sel]) )
            self.root_hists[ihist].SetBinError(self.root_hists[ihist].GetNbinsX() + 1, np.sqrt(np.sum(np.square(weight[full_sel]) )))

            self.bin_sum.append(np.sum([ self.root_hists[ihist].GetBinContent(ibin + 1) for ibin in range(self.root_hists[ihist].GetNbinsX()) ]))
            self.bin_norm.append(np.sum([ self.root_hists[ihist].GetBinContent(ibin) for ibin in range(self.root_hists[ihist].GetNbinsX() + 2) ]))

        self.dim1_underflow = np.sum(weight[inf_weight_mask & filter_cut & (obs_arrays[self.dim1.np_var] < self.dim1.bin_edges[0])])
        self.dim1_overflow = np.sum(weight[inf_weight_mask & filter_cut & (obs_arrays[self.dim1.np_var] >= self.dim1.bin_edges[len(self.dim1.bin_edges) - 1])])
        self.total = np.sum(self.bin_sum)
        self.norm = np.sum(weight[inf_weight_mask & filter_cut])
        return

    def flatten_hist(self):
        BinEdges = []
        self.flat_xAxisLabel = []
        self.flat_xAxisLabel_up = []
        self.flat_xAxisLabel_low = []
        start = 0.0
        for ihist, Hist in enumerate(self.root_hists):
            self.flat_xAxisLabel_low.append(start + Hist.GetXaxis().GetBinLowEdge(1))
            self.flat_xAxisLabel_up.append(start + Hist.GetXaxis().GetBinUpEdge(Hist.GetNbinsX()))
            self.flat_xAxisLabel.append(ROOT.TF1(f'labels{ihist}', 'x', Hist.GetXaxis().GetBinLowEdge(1), Hist.GetXaxis().GetBinUpEdge(Hist.GetNbinsX())))
            for ibin in range(Hist.GetNbinsX()):
                BinEdges.append(Hist.GetXaxis().GetBinLowEdge(ibin + 1) + start)

            start += Hist.GetXaxis().GetBinUpEdge(Hist.GetNbinsX()) - Hist.GetXaxis().GetBinLowEdge(1)

        BinEdges.append(start + Hist.GetXaxis().GetBinLowEdge(1))
        BinEdges = np.array(BinEdges)
        self.root_flat_hist = ROOT.TH1F(self.name + '_Merge_' + self.tag, self.name + '_Merge_' + self.tag, len(BinEdges) - 1, BinEdges)
        ibin_merge = 0
        for Hist in self.root_hists:
            for ibin in range(Hist.GetNbinsX()):
                self.root_flat_hist.SetBinContent(ibin_merge + 1, Hist.GetBinContent(ibin + 1))
                self.root_flat_hist.SetBinError(ibin_merge + 1, Hist.GetBinError(ibin + 1))
                ibin_merge += 1

    def fill_2Dhist(self, event_data, from_root, hist_name=None, weightarray=None, genWeight=''):
        if from_root:
            self.fill_2Dhist_from_root( event_data, hist_name=hist_name, genWeight=genWeight )
        else:
            self.fill_2Dhist_from_npz( event_data, hist_name=hist_name, weightarray=weightarray, genWeight=genWeight)

    def fill_2Dhist_from_root(self, tree, hist_name=None, genWeight=''):
        if hist_name is None:
            hist_name = self.name
        self.root_2Dhist = ROOT.TH2F(hist_name, hist_name, len(self.dim1.bin_edges) - 1, np.array(self.dim1.bin_edges), len(self.dim2.bin_edges[0]) - 1, np.array(self.dim2.bin_edges[0]))
        string_weight = ''
        if genWeight != '':
            string_weight = '*' + genWeight
        tree.Draw(f'{self.dim2.root_var}:{self.dim1.root_var}>>{hist_name}', self.root_cut + string_weight, 'colzgoff')
        self.bin_sum = [ np.sum([ self.root_2Dhist.GetBinContent(ibin1, ibin2) for ibin2 in range(1, self.root_2Dhist.GetNbinsY() + 1) ]) for ibin1 in range(1, self.root_2Dhist.GetNbinsX() + 1) ]
        self.bin_norm = [ np.sum([ self.root_2Dhist.GetBinContent(ibin1, ibin2) for ibin2 in range(self.root_2Dhist.GetNbinsY() + 2) ]) for ibin1 in range(1, self.root_2Dhist.GetNbinsX() + 1) ]
        self.dim1_underflow = np.sum([ self.root_2Dhist.GetBinContent(0, ibin2) for ibin2 in range(self.root_2Dhist.GetNbinsY() + 2) ])
        self.dim1_overflow = np.sum([ self.root_2Dhist.GetBinContent(self.root_2Dhist.GetNbinsX() + 1, ibin2) for ibin2 in range(self.root_2Dhist.GetNbinsY() + 2) ])
        self.total = np.sum(self.bin_sum)
        self.norm = np.sum([ np.sum([ self.root_2Dhist.GetBinContent(ibin1, ibin2) for ibin2 in range(self.root_2Dhist.GetNbinsY() + 2) ]) for ibin1 in range(self.root_2Dhist.GetNbinsX() + 2) ])
        return

    def fill_2Dhist_from_npz(self, files, hist_name=None, weightarray=None, genWeight=''):
        if hist_name is None:
            hist_name = self.name
        self.root_2Dhist = ROOT.TH2F(hist_name, hist_name, len(self.dim1.bin_edges) - 1, np.array(self.dim1.bin_edges), len(self.dim2.bin_edges[0]) - 1, np.array(self.dim2.bin_edges[0]))
        if isinstance(files, list):
            f_list = [ np.load(str(file_one), allow_pickle=True) for file_one in files ]
            obs_arrays = {}
            for key in list(f_list[0].keys()):
                if key != 'tracks' and key != 'charged':
                    obs_arrays[key] = np.concatenate([ f[key] for f in f_list ], axis=0)

        else:
            obs_arrays = np.load(files)
        if weightarray is None:
            weightarray = np.ones(len(obs_arrays['reco_ntrk']))
        if genWeight != '':
            if genWeight in obs_arrays.keys():
                weight = weightarray * obs_arrays[genWeight]
            else:
                print(f'{genWeight} is not a key in {files} skip it')
                weight = weightarray
        else:
            weight = weightarray
        filter_cut = filter_np_cut(obs_arrays, self.npy_cut)
        base_cut = (~(np.isinf(weight)) & filter_cut )
        for ihist in range(self.root_2Dhist.GetNbinsX()):
            ihist_cut = (obs_arrays[self.dim1.np_var] >= self.dim1.bin_edges[ihist]) & (obs_arrays[self.dim1.np_var] < self.dim1.bin_edges[ihist + 1]) 
            for ibin in range(self.root_2Dhist.GetNbinsY()):
                wgts = weight[base_cut & ihist_cut & (obs_arrays[self.dim2.np_var] >= self.dim2.bin_edges[ihist][ibin]) & (obs_arrays[self.dim2.np_var] < self.dim2.bin_edges[ihist][ibin + 1])]
                self.root_2Dhist.SetBinContent(ihist + 1, ibin + 1, np.sum(wgts) )
                self.root_2Dhist.SetBinError(ihist + 1, ibin + 1, np.sqrt(np.sum(np.square(wgts))) )

            wgts = weight[ base_cut &  ihist_cut & (obs_arrays[self.dim2.np_var] < self.dim2.bin_edges[ihist][0])]
            self.root_2Dhist.SetBinContent(ihist + 1, 0, np.sum(wgts) )
            self.root_2Dhist.SetBinError(ihist + 1, 0, np.sqrt(np.sum(np.square(wgts))) )

            wgts = weight[ base_cut & ihist_cut & (obs_arrays[self.dim2.np_var] >= self.dim2.bin_edges[ihist][len(self.dim2.bin_edges[ihist]) - 1])]
            self.root_2Dhist.SetBinContent(ihist + 1, self.root_2Dhist.GetNbinsY() + 1, np.sum(wgts) )
            self.root_2Dhist.SetBinError(ihist + 1, self.root_2Dhist.GetNbinsY() + 1, np.sqrt(np.sum(np.square(wgts))) )

            self.bin_sum.append(np.sum([ self.root_2Dhist.GetBinContent(ihist + 1, ibin + 1) for ibin in range(self.root_2Dhist.GetNbinsX()) ]))
            self.bin_norm.append(np.sum([ self.root_2Dhist.GetBinContent(ihist + 1, ibin) for ibin in range(self.root_2Dhist.GetNbinsX() + 2) ]))

        for ibin in range(self.root_2Dhist.GetNbinsY()):
            wgts = weight[ base_cut & (obs_arrays[self.dim1.np_var] < self.dim1.bin_edges[0]) & (obs_arrays[self.dim2.np_var] >= self.dim2.bin_edges[0][ibin]) & (obs_arrays[self.dim2.np_var] < self.dim2.bin_edges[0][ibin + 1])]
            self.root_2Dhist.SetBinContent(0, ibin + 1, np.sum(wgts) )
            self.root_2Dhist.SetBinError(0, ibin + 1, np.sqrt(np.sum(np.square(wgts)) ) )

            wgts = weight[ base_cut & (obs_arrays[self.dim1.np_var] >= self.dim1.bin_edges[len(self.dim1.bin_edges) - 1]) & (obs_arrays[self.dim2.np_var] >= self.dim2.bin_edges[0][ibin]) & (obs_arrays[self.dim2.np_var] < self.dim2.bin_edges[0][ibin + 1])] 
            self.root_2Dhist.SetBinContent(self.root_2Dhist.GetNbinsX() + 1, ibin + 1, np.sum(wgts) )
            self.root_2Dhist.SetBinError(self.root_2Dhist.GetNbinsX() + 1, ibin + 1, np.sqrt(np.sum(np.square(wgts))) )

        wgts = weight[ base_cut & (obs_arrays[self.dim1.np_var] < self.dim1.bin_edges[0]) & (obs_arrays[self.dim2.np_var] < self.dim2.bin_edges[0][0])]
        self.root_2Dhist.SetBinContent(0, 0, np.sum(wgts) )
        self.root_2Dhist.SetBinError(0, 0, np.sqrt(np.sum(np.square(wgts) ) ) )

        wgts = weight[base_cut & (obs_arrays[self.dim1.np_var] < self.dim1.bin_edges[0]) & (obs_arrays[self.dim2.np_var] >= self.dim2.bin_edges[0][len(self.dim2.bin_edges[0]) - 1])]
        self.root_2Dhist.SetBinContent(0, self.root_2Dhist.GetNbinsY() + 1, np.sum(wgts) )
        self.root_2Dhist.SetBinError(0, self.root_2Dhist.GetNbinsY() + 1, np.sqrt(np.sum(np.square(wgts))) )

        wgts = weight[base_cut & (obs_arrays[self.dim1.np_var] >= self.dim1.bin_edges[len(self.dim1.bin_edges) - 1]) & (obs_arrays[self.dim2.np_var] < self.dim2.bin_edges[0][0])]
        self.root_2Dhist.SetBinContent(self.root_2Dhist.GetNbinsX() + 1, 0, np.sum(wgts) )
        self.root_2Dhist.SetBinError(self.root_2Dhist.GetNbinsX() + 1, 0, np.sqrt(np.sum(np.square(wgts))) )

        wgts = weight[base_cut & (obs_arrays[self.dim1.np_var] >= self.dim1.bin_edges[len(self.dim1.bin_edges) - 1]) & (obs_arrays[self.dim2.np_var] >= self.dim2.bin_edges[0][len(self.dim2.bin_edges[0]) - 1])]
        
        self.root_2Dhist.SetBinContent(self.root_2Dhist.GetNbinsX() + 1, self.root_2Dhist.GetNbinsY() + 1, np.sum(wgts) )
        self.root_2Dhist.SetBinError(self.root_2Dhist.GetNbinsX() + 1, self.root_2Dhist.GetNbinsY() + 1, np.sqrt(np.sum(np.square(wgts))) )

        self.dim1_underflow = np.sum(weight[base_cut & (obs_arrays[self.dim1.np_var] < self.dim1.bin_edges[0])])
        self.dim1_overflow = np.sum(weight[base_cut & (obs_arrays[self.dim1.np_var] >= self.dim1.bin_edges[len(self.dim1.bin_edges) - 1])])
        self.total = np.sum(self.bin_sum)
        self.norm = np.sum(weight[(np.isinf(weight) != True) & filter_cut])
        return

    def get_hist_from_division(self, hist_list_1, hist_list_2):
        if self.root_hists is None:
            self.root_hists = [ ROOT.TH1F(hist_name, hist_name, len(self.dim2.bin_edges[ibin1]) - 1, np.array(self.dim2.bin_edges[ibin1])) for ibin1, hist_name in enumerate(self.root_hists_name) ]
        for ihist in range(len(self.dim1.bin_edges) - 1):
            for ibin in range(self.root_hists[ihist].GetNbinsX() + 2):
                bin_value_hist1 = hist_list_1.root_hists[ihist].GetBinContent(ibin)
                bin_error_hist1 = hist_list_1.root_hists[ihist].GetBinError(ibin)
                bin_value_hist2 = hist_list_2.root_hists[ihist].GetBinContent(ibin)
                bin_error_hist2 = hist_list_2.root_hists[ihist].GetBinError(ibin)
                if bin_value_hist2 != 0:
                    self.root_hists[ihist].SetBinContent(ibin, bin_value_hist1 / bin_value_hist2)
                    self.root_hists[ihist].SetBinError(ibin, np.sqrt(bin_error_hist1 ** 2 / bin_value_hist2 ** 2 + bin_value_hist1 ** 2 * bin_error_hist2 ** 2 / bin_value_hist2 ** 4))

            self.bin_sum.append(np.sum([ self.root_hists[ihist].GetBinContent(ibin + 1) for ibin in range(self.root_hists[ihist].GetNbinsX()) ]))
            self.bin_norm.append(np.sum([ self.root_hists[ihist].GetBinContent(ibin) for ibin in range(self.root_hists[ihist].GetNbinsX() + 2) ]))

        self.dim1_underflow = hist_list_1.dim1_underflow / hist_list_2.dim1_underflow if hist_list_2.dim1_underflow > 0 else 0
        self.dim1_overflow = hist_list_1.dim1_overflow / hist_list_2.dim1_overflow if hist_list_2.dim1_overflow > 0 else 0
        self.total = np.sum(self.bin_sum)
        self.norm = np.sum(self.bin_norm)
        return

    def get_hist_from_multiplication(self, hist1, hist2):
        if self.root_hists is None:
            self.root_hists = [ ROOT.TH1F(hist_name, hist_name, len(self.dim2.bin_edges) - 1, np.array(self.dim2.bin_edges[ibin1])) for ibin1, hist_name in enumerate(self.root_hists_name) ]
        for ihist in range(len(self.dim1.bin_edges) - 1):
            for ibin in range(self.root_hists[ihist].GetNbinsX() + 2):
                bin_value_hist1 = hist_list_1.root_hists[ihist].GetBinContent(ibin)
                bin_error_hist1 = hist_list_1.root_hists[ihist].GetBinError(ibin)
                bin_value_hist2 = hist_list_2.root_hists[ihist].GetBinContent(ibin)
                bin_error_hist2 = hist_list_2.root_hists[ihist].GetBinError(ibin)
                self.root_hists[ihist].SetBinContent(ibin, bin_value_hist1 * bin_value_hist2)
                self.root_hists[ihist].SetBinError(ibin, np.sqrt(bin_error_hist1 ** 2 * bin_value_hist2 ** 2 + bin_error_hist2 ** 2 * bin_value_hist1 ** 2))

            self.bin_sum.append(np.sum([ self.root_hists[ihist].GetBinContent(ibin + 1) for ibin in range(self.root_hists[ihist].GetNbinsX()) ]))
            self.bin_norm.append(np.sum([ self.root_hists[ihist].GetBinContent(ibin) for ibin in range(self.root_hists[ihist].GetNbinsX() + 2) ]))

        self.dim1_underflow = hist_list_1.dim1_underflow * hist_list_2.dim1_underflow
        self.dim1_overflow = hist_list_1.dim1_overflow * hist_list_2.dim1_overflow
        self.total = np.sum(self.bin_sum)
        self.norm = np.sum(self.bin_norm)
        return

    def read_hist_from_file(self, f):
        for hist_name in self.root_hists_name:
          if hist_name not in f.GetListOfKeys():
            print(hist_name,"is not found")
            return -1
        self.root_hists = [ f.Get(hist_name).Clone(hist_name + '_Clone') for hist_name in self.root_hists_name ]
        for ihist in range(len(self.root_hists)):
            self.bin_sum.append(np.sum([ self.root_hists[ihist].GetBinContent(ibin + 1) for ibin in range(self.root_hists[ihist].GetNbinsX()) ]))
            self.bin_norm.append(np.sum([ self.root_hists[ihist].GetBinContent(ibin) for ibin in range(self.root_hists[ihist].GetNbinsX() + 2) ]))

        self.total = np.sum(self.bin_sum)
        self.norm = np.sum(self.total)
        return 0

    def binwise_multiply(self, scalar_list=[], scalar_underflow=0.0, scalar_overflow=0.0):
        if len(scalar_list) != len(self.root_hists):
            raise ValueError(f'The length of the given list is not equal to the length of the histogram, {len(self.root_hists)}')
        for ihist, hist in enumerate(self.root_hists):
            for ibin in range(self.root_hists[ihist].GetNbinsX() + 2):
                hist.SetBinContent(ibin, hist.GetBinContent(ibin) * scalar_list[ihist])
                hist.SetBinError(ibin, hist.GetBinError(ibin) * scalar_list[ihist])

            self.bin_sum[ihist] *= scalar_list[ihist]
            self.bin_norm[ihist] *= scalar_list[ihist]

        self.total = np.sum(self.bin_sum)
        self.dim1_underflow *= scalar_underflow
        self.dim1_overflow *= scalar_overflow
        self.norm = self.total + self.dim1_underflow + self.dim1_overflow

    def divide_by_bin_width(self):
        for hist in self.root_hists:
            for ibin in range(hist.GetNbinsX()):
                hist.SetBinContent(ibin + 1, hist.GetBinContent(ibin + 1) / hist.GetXaxis().GetBinWidth(ibin + 1))
                hist.SetBinError(ibin + 1, hist.GetBinError(ibin + 1) / hist.GetXaxis().GetBinWidth(ibin + 1))

    def multiply(self, scalar):
        for ihist, hist in enumerate(self.root_hists):
            for ibin in range(hist.GetNbinsX() + 2):
                hist.SetBinContent(ibin, hist.GetBinContent(ibin) * scalar)
                hist.SetBinError(ibin, hist.GetBinError(ibin) * scalar)

            self.bin_sum[ihist] *= scalar
            self.bin_norm[ihist] *= scalar

        self.total *= scalar
        self.dim1_underflow *= scalar
        self.dim1_overflow *= scalar
        self.norm *= scalar

    def binwise_normalize(self):
        self.binwise_multiply(scalar_list=[ 1.0 / x for x in self.bin_norm ], scalar_underflow=1.0 / self.dim1_underflow if self.dim1_underflow > 0 else 0.0, scalar_overflow=1.0 / self.dim1_overflow if self.dim1_overflow > 0 else 0.0)

    def normalize(self):
        self.binwise_multiply(scalar_list=[[1.0 / self.norm]] * len(self.bin_norm), scalar_underflow=1.0 / self.norm, scalar_overflow=1.0 / self.norm)

    def binwise_multiply_2Dhist(self, scalar_list=[], scalar_underflow=0.0, scalar_overflow=0.0):
        if len(scalar_list) != self.root_2Dhist.GetNbinsX():
            raise ValueError(f'The length of the given list is not equal to the length of the histogram, {self.root_2Dhist.GetNbinsX()}')
        for ibin1 in range(1, self.root_2Dhist.GetNbinsX() + 1):
            for ibin2 in range(self.root_2Dhist.GetNbinsY() + 2):
                self.root_2Dhist.SetBinContent(ibin1, ibin2, self.root_2Dhist.GetBinContent(ibin1, ibin2) * scalar_list[ibin1 - 1])
                self.root_2Dhist.SetBinError(ibin1, ibin2, self.root_2Dhist.GetBinError(ibin1, ibin2) * scalar_list[ibin1 - 1])

            self.bin_sum[ibin1 - 1] *= scalar_list[ibin1 - 1]
            self.bin_norm[ibin1 - 1] *= scalar_list[ibin1 - 1]

        for ibin2 in range(self.root_2Dhist.GetNbinsY() + 2):
            self.root_2Dhist.SetBinContent(0, ibin2, self.root_2Dhist.GetBinContent(0, ibin2) * scalar_underflow)
            self.root_2Dhist.SetBinError(0, ibin2, self.root_2Dhist.GetBinError(0, ibin2) * scalar_underflow)
            self.root_2Dhist.SetBinContent(self.root_2Dhist.GetNbinsX() + 1, ibin2, self.root_2Dhist.GetBinContent(self.root_2Dhist.GetNbinsX() + 1, ibin2) * scalar_overflow)
            self.root_2Dhist.SetBinError(self.root_2Dhist.GetNbinsX() + 1, ibin2, self.root_2Dhist.GetBinError(self.root_2Dhist.GetNbinsX() + 1, ibin2) * scalar_overflow)

        self.total = np.sum(self.bin_sum)
        self.dim1_underflow *= scalar_underflow
        self.dim1_overflow *= scalar_overflow
        self.norm = self.total + self.dim1_underflow + self.dim1_overflow

    def binwise_normalize_2Dhist(self):
        self.binwise_multiply_2Dhist(scalar_list=[ 1.0 / x for x in self.bin_norm ], scalar_underflow=1.0 / self.dim1_underflow if self.dim1_underflow > 0 else 0.0, scalar_overflow=1.0 / self.dim1_overflow if self.dim1_overflow > 0 else 0.0)

    def normalize_2Dhist(self):
        self.binwise_multiply_2Dhist(scalar_list=[[1.0 / self.norm]] * len(self.bin_norm), scalar_underflow=1.0 / self.norm, scalar_overflow=1.0 / self.norm)

    def write_hist_list(self):
        for hist in self.root_hists:
            hist.Write()

    def write_2Dhist(self):
        self.root_2Dhist.Write()

    def write_all(self):
        self.write_2Dhist()
        self.write_hist_list()
