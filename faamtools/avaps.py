# -*- coding: utf-8 -*-
"""
Functions to process AVAPS dropsondes data
"""
import datetime
import netCDF4 as nc
import numpy as np

from . import DsFld, ObsData
from . import utils


def read_avaps_nc(fname, flds=None, time2datetime=False):
    """Read AVAPS dropsonde data from a NetCDF file.

    Open a NetCDF file and write data into `ObsData` instance of `DsFld` objects.
    Perform filtering of raw dropsonde data using `var_utils.filt_miss_row`

    Args:
    -----
        fname: str, file name
    Kwargs:
    -------
        flds: dict, names of variables to read from a dropsonde data file
              The default value is
              dict(time='time',hgt='alt',lon='lon',lat='lat',
                   u='u_wind',v='v_wind',wspd='wspd',wdir='wdir',
                   pres='pres',tdry='tdry',thta='theta',dhgt='dz',
                   tdew='dp',relh='rh',mixr='mr',thte='theta_e',thtv='theta_v')
        time2datetime: boolean, optional.
                       If True and `flds` dictionary contains 'time' key, convert array of
                       time values to `datetime.datetime` objects.
                       Requires `var_utils.timestr2datetime()` to parse time units.
                       Defaults to False.
    Returns:
    --------
        data: `ObsData` instance

    """

    if flds == None:
        flds = dict(time='time',hgt='alt',lon='lon',lat='lat',\
                    u='u_wind',v='v_wind',wspd='wspd',wdir='wdir',\
                    pres='pres',tdry='tdry',thta='theta',dhgt='dz',\
                    tdew='dp',relh='rh',mixr='mr',thte='theta_e',thtv='theta_v')

    f = nc.Dataset(fname)

    dum = ObsData()
    for i in flds:
        ncfld = f.variables[flds[i]]
        dum(**{i:DsFld(raw=ncfld[:],units=ncfld.units,long_name=ncfld.long_name)})

    flds_list = [ii for ii in flds] # to keep the order
    fil_list = utils.filt_miss_row(*[getattr(dum,ii).raw for ii in flds_list])

    data = ObsData()
    for i, j in enumerate(fil_list):
        data(**{flds_list[i]:DsFld(raw=getattr(dum,flds_list[i]).raw,\
                                   fil=j,\
                                   units=getattr(dum,flds_list[i]).units,\
                                   long_name=getattr(dum,flds_list[i]).long_name)})

    if time2datetime and 'time' in flds:
        if hasattr(data.time, 'units'):
            tbase, tstep_sec = utils.timestr2datetime(data.time.units)
            arr_sec2datetime = np.vectorize(lambda x: tbase + datetime.timedelta(seconds=x*tstep_sec))
            data.time.fil = arr_sec2datetime(data.time.fil)

    return data
