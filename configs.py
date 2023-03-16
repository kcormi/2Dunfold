from dataclasses import dataclass, field, InitVar
from typing import Optional, Union
import numpy as np
import yaml
import json

from unfold_utils import HistList

@dataclass
class ConfigBase:
    '''An Base Class for configuration objects'''

    @classmethod
    def from_json( cls, fname, keys=None):
        '''Return a class instance directly from a json file.'''
        return cls._from_file( fname, keys, use_json=True)

    @classmethod
    def from_yaml( cls, fname, keys=None):
        '''Return a class instance directly from a yaml file.'''
        return cls._from_file( fname, keys)

    @classmethod
    def _from_file( cls, fname, keys=None, use_json=False ):
        '''Return a class instance directly from a yaml or json file.

        :param fname: the name of the yaml file.
        :param keys: list of keys to be traversed to get to the object.

        if the class object is the only one in the yaml or json file, no keys need to be passed
        but if many objects are stored in the yaml file and you only want one, you can specify
        which one with an ordered list of keys to be traversed.
        '''
        load_func = yaml.safe_load if not use_json else json.load
        with open( fname, 'r') as fhandle:
            dct = load_func( fhandle)

        keys = [] if keys is None else keys
        for k in keys:
            dct = dct[k]
        print( dct )
        return cls( **dct )

@dataclass
class SampleTypePlotConfig(ConfigBase):
    '''A class which keeps track of plot settings for different result types
    To ensure consist styling for the same result type across different plots.'''
    color: int
    legend: str
    tag: str
    sys_reweight: bool = field(init=False)
    color_unfold: int

    @property
    def legend_refold(self):
        return f'{self.legend} refold' if not self.sys_reweight else f'{self.legend} reweight'

    @property
    def legend_unfold(self):
        return f'{self.legend} unfold' if not self.sys_reweight else f'{self.legend} reweight'

@dataclass
class HistConfig(ConfigBase):
    '''A class which keepts track of plotting configurations for a given histogram'''
    hist: HistList
    stat: int
    color: int
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
        return ('reco' in self.name) or ('refold' in self.name)

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

@dataclass
class BinningConfig(ConfigBase):
    '''A class for binning configuration'''
    min: InitVar[float] = None
    max: InitVar[float] = None
    nbins: InitVar[int] = None
    edges: list[float] = field(default_factory=list)

    def __post_init__(self, min, max, nbins):
         if len(self.edges) == 0:
             self.edges = np.arange( min, max, (max-min)/nbins )

    @property
    def nbins(self):
      return len(self.edges[:-1])

    @property
    def min(self):
      return self.edges[0]

    @property
    def max(self):
      return self.edges[-1]

@dataclass
class VarConfig(ConfigBase):
    '''A class which keeps track of variable information.'''
    name: str
    np_column: str = None
    root_branch: str = None
    shortname: str = None

    def __post_init__(self):
        if self.np_column is None:
            if self.root_branch is None:
                self.np_column = name
            else:
                self.np_column = self.root_branch
        if self.root_branch is None:
            self.root_branch = self.np_column
        if self.shortname is None:
            self.shortname = self.name

    @property
    def is_binned(self):
        return False

@dataclass
class BinnedVarConfig(VarConfig):
    '''A class which keeps track of variable configuration with a specific binning'''
    binning: Union[BinningConfig,dict] = field(default_factory=dict)

    def __post_init__(self):
        if isinstance(self.binning, dict):
            self.binning = BinningConfig( **self.binning )

    @property
    def is_binned(self):
        return True

@dataclass
class ObsConfig(ConfigBase):
   '''A class which keeps track of observables configuration, which are variables paired at reco and gen level'''
   reco: Union[VarConfig,dict]
   gen: Union[VarConfig,dict]
   require_extra_file: int = 0

   def __post_init__(self):
        if isinstance( self.reco, dict):
            if "binning" in self.reco:
                self.reco = BinnedVarConfig( **self.reco )
            else:
                self.reco = VarConfig( **self.reco )
        if isinstance( self.gen, dict):
            if "binning" in self.gen:
                self.gen = BinnedVarConfig( **self.gen )
            else:
                self.get = VarConfig( **self.gen )

if __name__=='__main__':

    for obs in ['nparticle', 'mass', 'spherocity', 'transverse_spherocity', 'thrust', 'transverse_thrust', 'broadening', 'isotropy']:
        varCfg = ObsConfig.from_yaml( 'config/observables.yml', keys=[obs] )
        print(varCfg)
