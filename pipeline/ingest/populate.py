import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from pipeline import psth


settings = {'reserve_jobs': True, 'suppress_errors': True, 'display_progress': False}

psth.UnitPsth.populate(**settings)

psth.PeriodSelectivity.populate(**settings)

psth.UnitSelectivity.populate(**settings)
