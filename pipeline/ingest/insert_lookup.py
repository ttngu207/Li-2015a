import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import re
import scipy.io as sio
from tqdm import tqdm
import pathlib
import datajoint as dj

from pipeline import lab, experiment, ephys, virus
from pipeline import parse_date, dict_to_hash


# ================== DEFINE LOOK-UP ==================

# ==================== Probe =====================
# Probe - NeuroNexus Silicon Probe
probe = 'A4x8-5mm-100-200-177'
lab.Probe.insert1({'probe': probe,
                   'probe_type': 'nn_silicon_probe'})
lab.Probe.Electrode.insert({'probe': probe, 'electrode': x} for x in range(1, 33))

electrode_group = {'probe': probe, 'electrode_group': 0}
electrode_group_member = [{**electrode_group, 'electrode': chn} for chn in range(1, 33)]
electrode_config_name = 'silicon32'  #
electrode_config_hash = dict_to_hash(
    {**electrode_group, **{str(idx): k for idx, k in enumerate(electrode_group_member)}})
lab.ElectrodeConfig.insert1({'probe': probe,
                             'electrode_config_hash': electrode_config_hash,
                             'electrode_config_name': electrode_config_name})
lab.ElectrodeConfig.ElectrodeGroup.insert1({'electrode_config_name': electrode_config_name,
                                            **electrode_group})
lab.ElectrodeConfig.Electrode.insert({'electrode_config_name': electrode_config_name, **member}
                                     for member in electrode_group_member)

# ==================== Brain Location =====================
experiment.BrainLocation.insert1({'brain_location_name': 'left_m2',
                                  'brain_area': 'M2',
                                  'hemisphere': 'left',
                                  'skull_reference': 'Bregma'})
experiment.BrainLocation.insert1({'brain_location_name': 'right_m2',
                                  'brain_area': 'M2',
                                  'hemisphere': 'right',
                                  'skull_reference': 'Bregma'})
experiment.BrainLocation.insert1({'brain_location_name': 'both_m2',
                                  'brain_area': 'M2',
                                  'hemisphere': 'both',
                                  'skull_reference': 'Bregma'})
experiment.BrainLocation.insert1({'brain_location_name': 'left_alm',
                                  'brain_area': 'ALM',
                                  'hemisphere': 'left',
                                  'skull_reference': 'Bregma'})
experiment.BrainLocation.insert1({'brain_location_name': 'right_alm',
                                  'brain_area': 'ALM',
                                  'hemisphere': 'right',
                                  'skull_reference': 'Bregma'})
experiment.BrainLocation.insert1({'brain_location_name': 'both_alm',
                                  'brain_area': 'ALM',
                                  'hemisphere': 'both',
                                  'skull_reference': 'Bregma'})
experiment.BrainLocation.insert1({'brain_location_name': 'left_pons',
                                  'brain_area': 'PONS',
                                  'hemisphere': 'left',
                                  'skull_reference': 'Bregma'})
experiment.BrainLocation.insert1({'brain_location_name': 'right_pons',
                                  'brain_area': 'PONS',
                                  'hemisphere': 'right',
                                  'skull_reference': 'Bregma'})
experiment.BrainLocation.insert1({'brain_location_name': 'both_pons',
                                  'brain_area': 'PONS',
                                  'hemisphere': 'both',
                                  'skull_reference': 'Bregma'})