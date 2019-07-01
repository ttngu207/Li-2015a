import datajoint as dj
import hashlib


def get_schema_name(name):
    try:
        return dj.config['custom']['{}.database'.format(name)]
    except KeyError:
        if name.startswith('ingest'):
            prefix = '{}_ingest_'.format(dj.config.get('database.user', 'map'))
        else:
            prefix = 'map_v1_'

    return prefix + name