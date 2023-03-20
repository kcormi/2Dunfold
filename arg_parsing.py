import yaml
from argparse import ArgumentTypeError

def get_all_observable_names( config ):
    with open(config["varunfold"], 'r') as f:
        obs_cfg = yaml.safe_load(f)
        all_keys = [ k for k in obs_cfg.keys() if k != 'binning' ]
    return all_keys

def check_observables_in_list( all_obs, obs_list ):
    for obs in obs_list:
        if not obs in obs_list:
            raise ValueError(f'The observable {obs} was not found in the list of observables in the configuration file. The list is {all_obs}')

def get_obs_binning_tuples( obs_list ):

    all_tuples = []
    for obs in obs_list:
        if ':' in obs:
           try: 
                obs, binning = map(str, obs.split(':') )
           except:
                raise ValueError(f'expected obs binning pair to be of the from <obs_name>:<binning_name> but found {obs}')
           all_tuples.append( tuple([obs, binning]) )
        else:
            all_tuples.append( tuple([obs, None ]) )
    return all_tuples
        

def parse_obs_args( config, obs_list):
    all_obs = get_all_observable_names(config)
    check_observables_in_list(all_obs, obs_list )
    obs_tuples = get_obs_binning_tuples( obs_list )
    return obs_tuples[0], obs_tuples[1]


def obs_pair( pair ):
    try:
        obs1, obs2 = map(str, pair.split(',') )
    except:
        raise ArgumentTypeError('obs_pairs must a pair of strings separated by a comma')

    return obs1, obs2
