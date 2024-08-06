from dataclasses import dataclass, field
from typing import Union

from configs import  ConfigBase
from unfold_utils import HistList
from plot_mplhep_utils import HistArray

@dataclass
class ResultPlotSettings(ConfigBase):
    '''A class which keeps track of plot settings for different result types'''
    color: Union[int,str]
    legend: str
    sys_reweight: bool = field(init=False)
    gen_reweight: bool = field(init=False)
    color_unfold: Union[int,str]

    @property
    def reweight_label(self):
        if self.gen_reweight:
            return "gen-reweight"
        elif self.sys_reweight:
            return "gen-reco-reweight"
        else:
            return ""


    @property
    def legend_refold(self):
        return f'{self.legend} refold' 

    @property
    def legend_unfold(self):
        return f'{self.legend} unfold'

    @property
    def tag(self):
        return self.legend

@dataclass
class HistConfig(ConfigBase):
    '''A class which keepts track of plotting configurations for a given histogram'''
    hist: Union[HistList,HistArray]
    stat: Union[int,HistList]
    color: Union[int,str]
    style: str
    legend: str

@dataclass
class PlotConfig(ConfigBase):
    '''A class which keeps track of plotting configuration for an entire plot, consisting of multiple histograms'''
    ref: HistConfig
    compare: list[HistConfig]
    name: str
    ratio: str

    @property
    def is_reco(self):
        return ('reco' in self.name) or ('refold' in self.name) or ('acc' in self.name)

    @property
    def hists(self):
        return [ c.hist for c in self.compare if c and c.hist]

    @property
    def colors(self):
        return [ c.color for c in self.compare if c and c.hist ]

    @property
    def legends(self):
        return [ c.legend for c in self.compare if c and c.hist ]

    @property
    def styles(self):
        return [ c.style for c in self.compare if c and c.hist ]
