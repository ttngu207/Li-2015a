
import numpy as np

import matplotlib.pyplot as plt

from pipeline.psth import TrialCondition
from pipeline.psth import UnitPsth
from pipeline import ephys, experiment


_plt_xmin = -3
_plt_xmax = 2


def _plot_spike_raster(ipsi, contra, vlines=[], ax=None, title=''):
    if not ax:
       fig, ax = plt.subplots(1, 1)

    ipsi_tr = ipsi['raster'][1]
    for i, tr in enumerate(set(ipsi['raster'][1])):
        ipsi_tr = np.where(ipsi['raster'][1] == tr, i, ipsi_tr)

    contra_tr = contra['raster'][1]
    for i, tr in enumerate(set(contra['raster'][1])):
        contra_tr = np.where(contra['raster'][1] == tr, i, contra_tr)

    ax.plot(ipsi['raster'][0], ipsi_tr, 'r.', markersize=1)
    ax.plot(contra['raster'][0], contra_tr + ipsi_tr.max() + 1, 'b.', markersize=1)

    for x in vlines:
        ax.axvline(x=x, linestyle='--', color='k')

    ax.set_axis_off()
    ax.set_xlim([_plt_xmin, _plt_xmax])
    ax.set_title(title)


def _plot_psth(ipsi, contra, vlines=[], ax=None, title=''):
    if not ax:
       fig, ax = plt.subplots(1, 1)

    ax.plot(contra['psth'][1], contra['psth'][0], 'b')
    ax.plot(ipsi['psth'][1], ipsi['psth'][0], 'r')

    for x in vlines:
        ax.axvline(x=x, linestyle='--', color='k')

    ax.set_ylabel('spikes/s')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_xlim([_plt_xmin, _plt_xmax])
    ax.set_xlabel('Time (s)')
    ax.set_title(title)


def plot_default_unit_psth(unit_key):
    """
    Default raster and PSTH plot for a specified unit - only {good, no early lick, correct trials} selected
    """

    hemi = (ephys.ProbeInsertion.InsertionLocation
            * experiment.BrainLocation & unit_key).fetch1('hemisphere')

    ipsi_hit_cond_key = (TrialCondition
                         & {'trial_condition_name': ('good_noearlylick_left_hit'
                                                     if hemi == 'left' else 'good_noearlylick_right_hit')}).fetch1('KEY')

    contra_hit_cond_key = (TrialCondition
                           & {'trial_condition_name': ('good_noearlylick_right_hit'
                                                       if hemi == 'left' else 'good_noearlylick_left_hit')}).fetch1('KEY')

    ipsi_hit_unit_psth = UnitPsth.get_plotting_data(
        unit_key, ipsi_hit_cond_key)

    contra_hit_unit_psth = UnitPsth.get_plotting_data(
        unit_key, contra_hit_cond_key)

    period_starts = (experiment.Period
                     & 'period in ("sample", "delay", "response")').fetch(
                         'period_start')

    fig, axs = plt.subplots(2, 1)

    _plot_spike_raster(ipsi_hit_unit_psth, contra_hit_unit_psth, ax=axs[0],
                       vlines=period_starts, title=f'Unit #: {unit_key["unit"]}')
    _plot_psth(ipsi_hit_unit_psth, contra_hit_unit_psth, vlines=period_starts, ax=axs[1])

