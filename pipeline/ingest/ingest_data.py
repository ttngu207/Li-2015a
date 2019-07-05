import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import re
import scipy.io as sio
from tqdm import tqdm
import pathlib
from decimal import Decimal
import datajoint as dj
import numpy as np

from pipeline import lab, experiment, ephys, virus
from pipeline import parse_date, dict_to_hash, time_unit_conversion_factor


# ==================== DEFINE CONSTANTS =====================

trial_type_str = ['HitR', 'HitL', 'ErrR', 'ErrL', 'NoLickR', 'NoLickL']
trial_type_mapper = {'HitR': ('hit', 'right'),
                     'HitL': ('hit', 'left'),
                     'ErrR': ('miss', 'right'),
                     'ErrL': ('miss', 'left'),
                     'NoLickR': ('ignore', 'right'),
                     'NoLickL': ('ignore', 'left')}

photostim_mapper = {1: 'PONS', 2: 'ALM'}

post_resp_tlim = 2  # a trial may last at most 2 seconds after response cue

task_protocol = {'task': 'audio delay', 'task_protocol': 1}

# ================== INGESTION OF DATA ==================

data_dir = pathlib.Path('data', 'data_structure')
data_files = data_dir.glob('*.mat')

for data_file in data_files:
    fname = data_file.stem
    subject_id = int(re.search('ANM\d+', fname).group().replace('ANM', ''))
    session_date = parse_date(re.search('_\d+', fname).group().replace('_', ''))

    session_key = (experiment.Session & {'subject_id': subject_id, 'session_date': session_date}).fetch1('KEY')

    print(f'-- Read {data_file} --')
    print(f'\tMatched: {session_key}')
    sess_data = sio.loadmat(data_file, struct_as_record = False, squeeze_me=True)['obj']

    # get time conversion factor - (-1) to take into account Matlab's 1-based indexing
    ts_time_conversion = time_unit_conversion_factor[
        sess_data.timeUnitNames[sess_data.timeSeriesArrayHash.value.timeUnit - 1]]
    trial_time_conversion = time_unit_conversion_factor[
        sess_data.timeUnitNames[sess_data.trialTimeUnit - 1]]

    # time-series data
    ts_tvec = sess_data.timeSeriesArrayHash.value.time * ts_time_conversion
    ts_trial = sess_data.timeSeriesArrayHash.value.trial
    lick_trace = sess_data.timeSeriesArrayHash.value.valueMatrix[:, 0]
    aom_input_trace = sess_data.timeSeriesArrayHash.value.valueMatrix[:, 1]
    laser_power = sess_data.timeSeriesArrayHash.value.valueMatrix[:, 2]

    # trial data
    trial_zip = zip(sess_data.trialIds, sess_data.trialStartTimes * trial_time_conversion,
                    sess_data.trialTypeMat[:6, :].T, sess_data.trialTypeMat[6, :].T,
                    sess_data.trialPropertiesHash.value[0] * trial_time_conversion,
                    sess_data.trialPropertiesHash.value[1] * trial_time_conversion,
                    sess_data.trialPropertiesHash.value[2] * trial_time_conversion,
                    sess_data.trialPropertiesHash.value[-1])

    photostims = (experiment.Photostim * experiment.BrainLocation & session_key)

    tr_start_idx = np.append(1, np.diff(ts_trial))
    photostim_start = [laser_power[ts_trial == tr_id] for tr_id in sess_data.trialIds]

    session_trials, behavior_trials, trial_events, photostim_trials, photostim_events = [], [], [], [], []
    for (tr_id, tr_start, trial_type_mtx, is_early_lick,
         sample_start, delay_start, response_start, photostim_type) in trial_zip:

        tkey = dict(session_key, trial=tr_id,
                    start_time=Decimal(tr_start),
                    stop_time=Decimal(tr_start + response_start + post_resp_tlim))
        session_trials.append(tkey)

        trial_type = np.array(trial_type_str)[trial_type_mtx.astype(bool)]
        if len(trial_type) == 1:
            outcome, trial_instruction = trial_type_mapper[trial_type[0]]
        else:
            outcome, trial_instruction = 'non-performing', 'non-performing'

        bkey = dict(tkey, **task_protocol,
                    trial_instruction=trial_instruction,
                    outcome=outcome,
                    early_lick='early' if is_early_lick else 'no early')
        behavior_trials.append(bkey)

        for etype, etime in zip(('sample', 'delay', 'go'), (sample_start, delay_start, response_start)):
            if not np.isnan(etime):
                trial_events.append(dict(tkey, trial_event_id=len(trial_events)+1,
                                         trial_event_type=etype, trial_event_time=etime))

        if photostim_type != 0:
            pkey = dict(tkey)
            photostim_trials.append(pkey)
            if photostim_type in (1, 2):
                photostim_key = (photostims & {'brain_area': photostim_mapper[photostim_type.astype(int)]}).fetch1('KEY')
                stim_power = laser_power[ts_trial == tr_id]
                photostim_events.append(dict(pkey, **photostim_key, photostim_event_id=len(photostim_events)+1, power=stim_power.max()))


    # ---- insert TRIALS ----
    experiment.SessionTrial.insert(session_trials, allow_direct_insert=True)
    experiment.BehaviorTrial.insert(behavior_trials, allow_direct_insert=True)
    experiment.PhotostimTrial.insert(photostim_trials, allow_direct_insert=True)
    experiment.TrialEvent.insert(trial_events, allow_direct_insert=True)
    experiment.PhotostimEvent.insert(photostim_events, allow_direct_insert=True)









