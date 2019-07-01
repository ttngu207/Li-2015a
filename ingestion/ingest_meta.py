import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import re
from datetime import datetime

import numpy as np
from decimal import Decimal
import scipy.io as sio
import pandas as pd
from tqdm import tqdm
import pathlib
import datajoint as dj

from . import get_schema_name


lab = dj.create_virtual_module('lab', get_schema_name('lab'))
experiment = dj.create_virtual_module('experiment', get_schema_name('experiment'))

# ================== Ingestion of Metadata ==================

meta_data_dir = pathlib.Path('data', 'meta_data')

meta_data_files = meta_data_dir.glob('*.mat')
for meta_data_file in meta_data_files:
    print(f'-- Read {meta_data_file} --')
    meta_data = sio.loadmat(meta_data_file, struct_as_record = False, squeeze_me=True)['meta_data']

    # ---- subject ----
    subject_key = dict(subject_id=int(re.search('\d+', meta_data.animalID).group()),
                       cage_number=-1,
                       date_of_birth=datetime.strptime(meta_data.dateOfBirth, '%Y%m%d'),
                       sex=meta_data.sex.upper(),
                       animal_source=meta_data.animalSource)

    # ---- subject gene modification ----















