"""module to format and process values from a yt dataset"""
import os
import numpy as np
from .hdf5yt import getFields


def get_plotTitles(ds, comment='', extras={}):
    """builds a pair or strings with information from the run.
    Namely the match radius and size, resolution of the grid,
    number of species, timestamp and profile used

    Args:
        ds(yt.dataset): loaded yt dataset.

    Returns:
        (str list): two axes titles for figure.

    """
    # match size and location radius
    mrad = ds.parameters['x_match']**2
    mrad += ds.parameters['y_match']**2
    mrad += ds.parameters['z_match']**2
    mrad = np.sqrt(mrad)/1e5  # factor to km
    msize = ds.parameters['r_match_outer']/1e5
    # number of species
    _, sps = getFields(ds.field_list)
    # figure out resolution only on 'x'
    res = ds.domain_width[0].value
    res /= ds.parameters['nxb']
    res /= ds.parameters['nblockx']
    maxres = res/2**(ds.parameters['lrefine_max']-1)/1e5
    norres = res/2**(ds.parameters['refnogenlevel']-1)/1e5
#     name = os.path.basename(ds.parameters['initialwdfile'][:-4])
    name = os.path.basename(os.path.dirname(ds.parameters['log_file']))
    datb = "{}\n".format(name)
    datb += "Max(nonENUC) res:{:2.0f}({:2.0f})km\n".format(maxres, norres)
    datb += "Match(radius):{:5.0f}({:2.0f})km\n".format(mrad, msize)
    datb += "Species: {:16d}".format(len(sps))
    # add time and max speed
    datb2 = comment + '\n'
    datb2 += "Time: {:1.7f} s\n".format(float(ds.current_time))
    if 'maxspeed' in extras:
        datb2 += "Max Speed: {:11.0f} km/s".format(extras['maxspeed']/1e5)
    else:
        for k, v in extras.items():
            datb2 += "{}: {:.3E}".format(k, v)
    return datb, datb2
