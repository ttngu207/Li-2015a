#!/usr/bin/env python3
import os

import sys
from datetime import datetime
from dateutil.tz import tzlocal
import pytz
import re
import numpy as np
import pandas as pd

from pipeline import (lab, experiment, ephys, psth, tracking)
import pynwb
from pynwb import NWBFile, NWBHDF5IO

# ============================== SET CONSTANTS ==========================================
nwb_output_dir = os.path.join('data', 'NWB 2.0')
session_time = datetime.strptime('00:00:00', '%H:%M:%S').time()  # no precise time available
hardware_filter = 'Bandpass filtered 300-6K Hz'


def export_to_nwb(session_key):

    this_session = (experiment.Session & session_key).fetch1()

    # ===============================================================================
    # ============================== META INFORMATION ===============================
    # ===============================================================================

    # -- NWB file - a NWB2.0 file for each session
    nwbfile = NWBFile(identifier='_'.join(
        [str(this_session['subject_id']),
         this_session['session_date'].strftime('%Y-%m-%d'),
         str(this_session['session'])]),
        session_description='',
        session_start_time=datetime.combine(this_session['session_date'], session_time),
        file_create_date=datetime.now(tzlocal()),
        experimenter=this_session['username'],
        institution='Janelia Research Campus')
    # -- subject
    subj = (lab.Subject & session_key).fetch1()
    nwbfile.subject = pynwb.file.Subject(
        subject_id=str(this_session['subject_id']),
        genotype=' x '.join((lab.Subject.GeneModification
                             & subj).fetch('gene_modification')),
        sex=subj['sex'],
        species=subj['species'],
        date_of_birth=subj['date_of_birth'])

    # ===============================================================================
    # ======================== EXTRACELLULAR & CLUSTERING ===========================
    # ===============================================================================

    """
    In the event of multiple probe recording (i.e. multiple probe insertions), the clustering results 
    (and the associated units) are associated with the corresponding probe. 
    Each probe insertion is associated with one ElectrodeConfiguration (which may define multiple electrode groups)
    """

    dj_insert_location = ephys.ProbeInsertion.InsertionLocation * experiment.BrainLocation

    for probe_insertion in ephys.ProbeInsertion & session_key:
        electrode_config = (lab.ElectrodeConfig & probe_insertion).fetch1()

        electrode_groups = {}
        for electrode_group in lab.ElectrodeConfig.ElectrodeGroup & electrode_config:
            electrode_groups[electrode_group['electrode_group']] = nwbfile.create_electrode_group(
                name=electrode_config['electrode_config_name'] + '_g' + str(electrode_group['electrode_group']),
                description = 'N/A',
                device = nwbfile.create_device(name=electrode_config['probe']),
                location = '; '.join([f'{k}: {str(v)}' for k, v in
                                      (dj_insert_location & session_key).fetch1().items()
                                      if k not in dj_insert_location.primary_key]))

        for chn in (lab.ElectrodeConfig.Electrode * lab.Probe.Electrode & electrode_config).fetch(as_dict=True):
            nwbfile.add_electrode(id=chn['electrode'],
                                  group=electrode_groups[chn['electrode_group']],
                                  filtering=hardware_filter,
                                  imp=-1.,
                                  x=chn['x_coord'] if chn['x_coord'] else np.nan,
                                  y=chn['y_coord'] if chn['y_coord'] else np.nan,
                                  z=chn['z_coord'] if chn['z_coord'] else np.nan,
                                  location=electrode_groups[chn['electrode_group']].location)

        # --- unit spike times ---
        nwbfile.add_unit_column(name='quality', description='unit quality from clustering')
        nwbfile.add_unit_column(name='posx', description='estimated x position of the unit relative to probe (0,0)')
        nwbfile.add_unit_column(name='posy', description='estimated y position of the unit relative to probe (0,0)')
        nwbfile.add_unit_column(name='amp', description='unit amplitude')
        nwbfile.add_unit_column(name='snr', description='unit signal-to-noise')
        nwbfile.add_unit_column(name='cell_type', description='cell type (e.g. fast spiking or pyramidal)')

        for unit in (ephys.Unit * ephys.UnitCellType & probe_insertion).fetch(as_dict=True):
            # make an electrode table region (which electrode(s) is this unit coming from)
            nwbfile.add_unit(id=unit['unit'],
                             electrodes=unit['electrode'],
                             electrode_group=electrode_groups[unit['electrode_group']],
                             quality = unit['unit_quality'],
                             posx=unit['unit_posx'],
                             posy=unit['unit_posy'],
                             amp = unit['unit_amp'],
                             snr = unit['unit_snr'],
                             cell_type=unit['cell_type'],
                             spike_times=unit['spike_times'],
                             waveform_mean=np.mean(unit['waveform'], axis=0),
                             waveform_sd=np.std(unit['waveform'], axis=0))

    # ===============================================================================
    # ============================= BEHAVIOR TRACKING ===============================
    # ===============================================================================
    tracking_data = ((tracking.LickTrace & session_key).fetch1()
                     if tracking.LickTrace & session_key else None)
    
    if tracking_data:
        behav_acq = pynwb.behavior.BehavioralTimeSeries(name = 'lick_trace')
        nwbfile.add_acquisition(behav_acq)
        [tracking_data.pop(k) for k in tracking.LickTrace.primary_key]
        timestamps = tracking_data.pop('lick_trace_timestamps')
        for b_k, b_v in tracking_data.items():
            behav_acq.create_timeseries(name = b_k,
                                        unit = 'a.u.',
                                        conversion = 1.0,
                                        data = b_v,
                                        timestamps=timestamps)

    # ===============================================================================
    # ============================= PHOTO-STIMULATION ===============================
    # ===============================================================================
    for photostim in experiment.Photostim * experiment.BrainLocation * lab.PhotostimDevice & session_key:

        stim_device = (nwbfile.get_device(photostim['photostim_device'])
                       if photostim['photostim_device'] in nwbfile.devices
                       else nwbfile.create_device(name=photostim['photostim_device']))

        stim_site = pynwb.ogen.OptogeneticStimulusSite(
            name=photostim['brain_location_name'],
            device=stim_device,
            excitation_lambda=float(photostim['excitation_wavelength']),
            location='; '.join([f'{k}: {str(v)}' for k, v in photostim.items()
                                if k in dj_insert_location.heading.names and k not in dj_insert_location.primary_key]),
            description='')
        nwbfile.add_ogen_site(stim_site)

    # ===============================================================================
    # =============================== BEHAVIOR TRIALS ===============================
    # ===============================================================================

    # =============== TrialSet ====================
    # NWB 'trial' (of type dynamic table) by default comes with three mandatory attributes:
    #                                                                       'id', 'start_time' and 'stop_time'.
    # Other trial-related information needs to be added in to the trial-table as additional columns (with column name
    # and column description)

    dj_trial = experiment.SessionTrial * experiment.BehaviorTrial
    skip_adding_columns = experiment.Session.primary_key + ['trial_uid']

    if experiment.SessionTrial & session_key:
        # Get trial descriptors from TrialSet.Trial and TrialStimInfo
        trial_columns = [{'name': tag,
                          'description': re.sub('\s+:|\s+', ' ', re.search(
                              f'(?<={tag})(.*)', str(dj_trial.heading)).group()).strip()}
                         for tag in dj_trial.heading.names
                         if tag not in skip_adding_columns + ['start_time', 'stop_time']]

        # Add new table columns to nwb trial-table for trial-label
        for c in trial_columns:
            nwbfile.add_trial_column(**c)

        # Add entry to the trial-table
        for trial in (dj_trial & session_key).fetch(as_dict=True):
            trial['start_time'] = float(trial['start_time'])
            trial['stop_time'] = float(trial['stop_time']) if trial['stop_time'] else 5.0
            [trial.pop(k) for k in skip_adding_columns]
            nwbfile.add_trial(**trial)

    # ===============================================================================
    # =============================== TRIAL-SEGMENTED DATA ==========================
    # ===============================================================================

    # =============== Write NWB 2.0 file ===============
    save_file_name = ''.join([nwbfile.identifier, '.nwb'])
    if not os.path.exists(nwb_output_dir):
        os.makedirs(nwb_output_dir)
    with NWBHDF5IO(os.path.join(nwb_output_dir, save_file_name), mode = 'w') as io:
        io.write(nwbfile)
        print(f'Write NWB 2.0 file: {save_file_name}')






