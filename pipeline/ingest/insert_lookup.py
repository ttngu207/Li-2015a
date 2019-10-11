import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
import itertools

from pipeline import lab, experiment, psth
from pipeline import dict_to_hash


# ================== DEFINE LOOK-UP ==================

# ==================== Project =====================

experiment.Project.insert([('li2015', '', 'doi:10.1038/nature14178'),
                           ('lidaie2016', '', 'doi:10.1038/nature17643')],
                          skip_duplicates=True)

# ==================== Probe =====================
# Probe - NeuroNexus Silicon Probe
probe = 'A4x8-5mm-100-200-177'
lab.Probe.insert1({'probe': probe,
                   'probe_type': 'nn_silicon_probe'}, skip_duplicates=True)
lab.Probe.Electrode.insert(({'probe': probe, 'electrode': x} for x in range(1, 33)), skip_duplicates=True)

electrode_group = {'probe': probe, 'electrode_group': 0}
electrode_group_member = [{**electrode_group, 'electrode': chn} for chn in range(1, 33)]
electrode_config_name = 'silicon32'  #
electrode_config_hash = dict_to_hash(
    {**electrode_group, **{str(idx): k for idx, k in enumerate(electrode_group_member)}})
lab.ElectrodeConfig.insert1({'probe': probe,
                             'electrode_config_hash': electrode_config_hash,
                             'electrode_config_name': electrode_config_name}, skip_duplicates=True)
lab.ElectrodeConfig.ElectrodeGroup.insert1({'electrode_config_name': electrode_config_name,
                                            **electrode_group}, skip_duplicates=True)
lab.ElectrodeConfig.Electrode.insert(({'electrode_config_name': electrode_config_name, **member}
                                     for member in electrode_group_member), skip_duplicates=True)

# ==================== Brain Location =====================
brain_locations = [{'brain_location_name': 'left_m2',
                    'brain_area': 'M2',
                    'hemisphere': 'left',
                    'skull_reference': 'Bregma'},
                   {'brain_location_name': 'right_m2',
                    'brain_area': 'M2',
                    'hemisphere': 'right',
                    'skull_reference': 'Bregma'},
                   {'brain_location_name': 'both_m2',
                    'brain_area': 'M2',
                    'hemisphere': 'both',
                    'skull_reference': 'Bregma'},
                   {'brain_location_name': 'left_alm',
                    'brain_area': 'ALM',
                    'hemisphere': 'left',
                    'skull_reference': 'Bregma'},
                   {'brain_location_name': 'right_alm',
                    'brain_area': 'ALM',
                    'hemisphere': 'right',
                    'skull_reference': 'Bregma'},
                   {'brain_location_name': 'both_alm',
                    'brain_area': 'ALM',
                    'hemisphere': 'both',
                    'skull_reference': 'Bregma'},
                   {'brain_location_name': 'left_pons',
                    'brain_area': 'PONS',
                    'hemisphere': 'left',
                    'skull_reference': 'Bregma'},
                   {'brain_location_name': 'right_pons',
                    'brain_area': 'PONS',
                    'hemisphere': 'right',
                    'skull_reference': 'Bregma'},
                   {'brain_location_name': 'both_pons',
                    'brain_area': 'PONS',
                    'hemisphere': 'both',
                    'skull_reference': 'Bregma'}]
experiment.BrainLocation.insert(brain_locations, skip_duplicates=True)

# ==================== Photostim Trial Condition =====================

stim_locs = ['left_alm', 'right_alm', 'both_alm']
stim_periods = [None, 'sample', 'early_delay', 'middle_delay']

trial_conditions = []
for loc in stim_locs:
    for instruction in (None, 'left', 'right'):
        for period, stim_dur in itertools.product(stim_periods, (0.5, 0.8)):
            condition = {'trial_condition_name': '_'.join(filter(None, ['all', 'noearlylick', loc,
                                                                        period, str(stim_dur), 'stim', instruction])),
                         'trial_condition_func': '_get_trials_include_stim',
                         'trial_condition_arg': {
                             **{'_outcome': 'ignore',
                                'task': 'audio delay',
                                'task_protocol': 1,
                                'early_lick': 'no early',
                                'brain_location_name': loc},
                             **({'trial_instruction': instruction} if instruction else {'_trial_instruction': 'non-performing'}),
                             **({'photostim_period': period, 'duration': stim_dur} if period else dict())}}
            trial_conditions.append(condition)

psth.TrialCondition.insert_trial_conditions(trial_conditions)

