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

def parse_obs_args( config, obs_list):
    all_obs = get_all_observable_names(config)
    check_observables_in_list(all_obs, obs_list )
    return obs_list[0], obs_list[1]


def obs_pair( pair ):
    try:
        obs1, obs2 = map(str, pair.split(',') )
    except:
        raise ArgumentTypeError('obs_pairs must a pair of strings separated by a comma')

    return obs1, obs2
