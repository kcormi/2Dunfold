from dataclasses import dataclass, field, InitVar
from typing import Optional, Union
import numpy as np
import yaml
import json

@dataclass
class ConfigBase:
    '''An Base Class for configuration objects'''

    def __post_init__(self):
        pass

    @classmethod
    def from_json( cls, fname, keys=None, **kwargs):
        '''Return a class instance directly from a json file.'''
        return cls._from_file( fname, keys, use_json=True, **kwargs)

    @classmethod
    def from_yaml( cls, fname, keys=None, **kwargs):
        '''Return a class instance directly from a yaml file.'''
        return cls._from_file( fname, keys, **kwargs)

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
class BinningConfig(ConfigBase):
    '''A class for binning configuration'''
    min: InitVar[float] = None
    max: InitVar[float] = None
    nbins: InitVar[int] = None
    edges: list[float] = field(default_factory=list)

    def __post_init__(self, min, max, nbins):
         super().__post_init__()
         if len(self.edges) == 0:
             self.edges = np.linspace( min, max, num=nbins+1 )

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
    np_var: str = None
    root_var: str = None
    shortname: str = None

    def __post_init__(self):
        super().__post_init__()
        if self.np_var is None:
            if self.root_var is None:
                self.np_var = name
            else:
                self.np_var = self.root_var
        if self.root_var is None:
            self.root_var = self.np_var
        if self.shortname is None:
            self.shortname = self.name

    @property
    def is_binned(self):
        return self.getattr('edges',None) is not None

@dataclass
class BinnedVarConfig(BinningConfig,VarConfig):
    '''A class which keeps track of variable configuration with a specific binning'''

    pass

@dataclass
class ObsConfig(ConfigBase):
    '''A class which keeps track of observables configuration, which are variables paired at reco and gen level'''
    reco: Union[VarConfig,dict]
    gen: Union[VarConfig,dict]
    require_extra_file: int = 0

    def __post_init__(self):
         super().__post_init__()
         if isinstance( self.reco, dict):
             if "edges" in self.reco or "nbins" in self.reco:
                 self.reco = BinnedVarConfig( **self.reco )
             else:
                 self.reco = VarConfig( **self.reco )
         if isinstance( self.gen, dict):
             if "edges" in self.gen or "nbins" in self.gen:
                 self.gen = BinnedVarConfig( **self.gen )
             else:
                 self.get = VarConfig( **self.gen )

    def __getitem__(self, key):
        if key == 'reco':
            return self.reco
        elif key == 'gen':
            return self.gen
        else:
            raise KeyError(f'No key {key} if ObsConfig, keys are "reco" and "gen"')

    @classmethod
    def _get_binning_from_dict( cls, binning, dct ):
        binning_meta_name = 'binnings'
        full_binning_dct = dct.get(binning_meta_name, None )
        if full_binning_dct is None:
            raise ValueError(f'could not find name binnings in the configuration, they should be under the top level heading {binning_meta_name}. But found only: {dct.keys()}')

        this_binning_dct = full_binning_dct.get( binning, None)
        if this_binning_dct is None:
            raise ValueError(f'could not find a binning by the name {binning} in the configuration, only found options {full_binning_dct.keys()}')
        return this_binning_dct

    @classmethod
    def _from_file(cls, fname, keys=None, use_json=False, binning=None):
        '''Load from file, but if an optional binning parameter is passed,
            find the appropriate binning object and use that instead of the default
            binnings.'''
        if binning is None:
            return super()._from_file( fname, keys=keys, use_json=use_json)
        else:
            load_func = yaml.safe_load if not use_json else json.load
            with open( fname, 'r') as fhandle:
                dct = load_func( fhandle)
            binning_dct = cls._get_binning_from_dict( binning, dct )


            keys = [] if keys is None else keys
            for k in keys:
                dct = dct[k]

            for bin_key in BinningConfig.__annotations__.keys():
                if bin_key in dct:
                    dct.pop(bin_key)

            for var_type in ['reco','gen']:
                dct[var_type].update( binning_dct[var_type] )

            return cls( **dct )


@dataclass
class HistDim(BinnedVarConfig):
    '''Class containing information about the dimension of a histogram,
        this may include bins of one variable nested inside another variable.'''
    inner_dim: BinnedVarConfig = None
    underflow: float = 0.
    overflow: float = 0.

    def __post_init__(self, *args):
        super().__post_init__(*args)
        inner_dim_factor = 0
        if self.inner_dim is not None:
            inner_dim_factor = max(1, self.inner_dim.nbins )
            self.edges = [self.edges] * inner_dim_factor

