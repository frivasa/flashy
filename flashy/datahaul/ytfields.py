"""custom fields for yt.
method() describe operations over the ds
_method fills the metadata for yt:
units, take_log, sampling_type
"""
from flashy.utils import np
_speed = {'units': 'cm/s', 'take_log': False, 'sampling_type': 'cell'}


def speed(field, data):
    vx = data['flash', 'velx']
    vy = data['flash', 'vely']
    vz = data['flash', 'velz']
    spd = np.sqrt(vx*vx + vy*vy, vz*vz)
    return spd
