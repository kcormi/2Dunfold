from dataclasses import dataclass, field
from typing import Union

from configs import  ConfigBase
from unfold_utils import HistList


@dataclass
class ResultPlotSettings(ConfigBase):
    '''A class which keeps track of plot settings for different result types'''
    color: Union[int,str]
    legend: str
    sys_reweight: bool = field(init=False)
    color_unfold: Union[int,str]

    @property
    def legend_refold(self):
        return f'{self.legend} refold' if not self.sys_reweight else f'{self.legend} reweight'

    @property
    def legend_unfold(self):
        return f'{self.legend} unfold' if not self.sys_reweight else f'{self.legend} reweight'

    @property
    def tag(self):
        return self.legend

@dataclass
class HistConfig(ConfigBase):
    '''A class which keepts track of plotting configurations for a given histogram'''
    hist: HistList
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
        return [ c.hist for c in self.compare if c.hist ]

    @property
    def colors(self):
        return [ c.color for c in self.compare if c.hist ]

    @property
    def legends(self):
        return [ c.legend for c in self.compare if c.hist ]

    @property
    def styles(self):
        return [ c.style for c in self.compare if c.hist ]
