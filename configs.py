from dataclasses import dataclass, field

@dataclass
class SampleTypePlotConfig:
    '''A class which keeps track of plot settings for different result types
    To ensure consist styling for the same result type across different plots'''
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
class HistConfig:
    '''A class which keepts track of plotting configurations for a given histogram'''
    hist: HistList
    stat: int
    color: int
    style: str
    legend: str

@dataclass
class PlotConfig:
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
def VarConfig:
   '''A class which keeps track of a variables configuration'''
   reco: str
   gen: str
   reco_name: str
   gen_name: str
   reco_shortname: str = field(init=False)
   gen_shortname: str = field(init=False)
   reco_key: str = field(init=False)
   gen_key: str = field(init=False)
   nbinsreco: int
   minreco: float
   maxreco: float
   binedgesreco: list[float] = field(default_factory=list)
   nbinsgen: int
   mingen: float
   maxgen: float
   binedgesgen: list[float] = field(default_factor=list)

  def __post_init__(self):
        if len(binedgesreco) == 0:
            self.binedgesreco = np.arange( self.minreco, self.maxreco, (self.maxreco-self.minreco)/self.nbinsreco )
        if len(self.binedgesgen) == 0:
            self.binedgesreco = np.arange( self.mingen, self.maxgen, (self.maxgen-self.mingen)/self.nbinsgen )
        if self.gen_shortname is None:
            self.gen_shortname = self.gen_name
        if self.reco_shortname is None:
            self.reco_shortname = self.reco_name

  @property
  def nbinsreco(self):
    return len(self.binedgesreco[:-1])

  @property
  def minreco(self):
    return self.binedgesreco[0]

  @property
  def maxreco(self):
    return self.binedgesreco[-1]

  @property
  def nbinsgen(self):
    return len(self.binedgesgen[:-1])

  @property
  def mingen(self):
    return self.binedgesgen[0]

  @property
  def maxgen(self):
    return self.binedgesgem[-1]

