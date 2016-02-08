# -*- coding: utf-8 -*-
"""
Functions to read FAAM core processed data
"""
import datetime
import netCDF4 as nc
import numpy as np

from . import FaamFld, ObsData
from . import utils


def read_core_nc(fname, flds=None, time2datetime=False, calc_wspd=True, calc_wdir=True):
    """Read core FAAM data from a NetCDF file.

    Open a NetCDF file and write data into `ObsData` instance of `FaamFld` objects

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
        calc_wspd: boolean, optional.
                   If True and `flds` dictionary contains 'u' and 'v' keys,
                   add calculate wind speed and add it to the `ObsData` instance.
                   Requires `var_utils.uv2wspd`. Defaults to True.
        calc_wdir: boolean, optional.
                   If True and `flds` dictionary contains 'u' and 'v' keys,
                   add calculate wind direction (degrees from North) and add it to the `ObsData` instance.
                   Requires `var_utils.uv2wdir`. Defaults to True.

    Returns:
    --------
        data: `ObsData` instance

    """

    if flds == None:
        flds = dict(time='Time',hgt='HGT_RADR',lon='LON_GPS',lat='LAT_GPS',ang='TRCK_GIN',\
                    uturb='U_C',vturb='V_C',wturb='W_C',\
                    u='U_NOTURB',v='V_NOTURB',\
                    pres='PS_RVSM',temp='TAT_DI_R')

    with nc.Dataset(fname) as f:

        data = ObsData()
        for i in flds:
            ncfld = f.variables[flds[i]]
            ncdata = ncfld[:]
            try:
                ncflag = f.variables[flds[i]+'_FLAG']
                ncdata[ncflag[:]!=0] = float('nan')
            except KeyError:
                pass

            data(**{i:FaamFld(val=ncdata[:].squeeze(),units=ncfld.units,long_name=ncfld.long_name)})

    if hasattr(data,'u') and hasattr(data,'v'):
        if calc_wspd:
            data.wspd = FaamFld(utils.uv2wspd(data.u.val,data.v.val),data.u.units,'wind speed derived from aircraft instruments and GIN')
        if calc_wdir:
            data.wdir = FaamFld(utils.uv2wdir(data.u.val,data.v.val),'deg','wind direction')

    if time2datetime and 'time' in flds:
        if hasattr(data.time, 'units'):
            tbase, tstep_sec = utils.timestr2datetime(data.time.units)
            arr_sec2datetime = np.vectorize(lambda x: tbase + datetime.timedelta(seconds=int(x)*tstep_sec))
            data.time.val = arr_sec2datetime(data.time.val)

    return data


def parse_profiles_runs_info(text_file_name, daystr='', timesorted=True):
    """Parse text file containing flight profiles and runs times
       a.k.a. FAAM sawtooth summary.

       Args:
       -----
           test_file_name: str, file name

       Kwargs:
       -------
           daystr: str, in format '%Y%m%d' date of observations, e.g. 20130326.
                   Defaults to an empty string.
           timesorted: bool, if True, sorts the tuples time-wise

       Returns:
       --------
           list of tuples like (name, start_time, finish_time)

       Example of sawtooth summary text file:
       --------------------------------------
           Profile 1
           111420
           111822
           Profile 2
           121459
           123830
           Profile 3
           131540
           131654
    """
    with open(text_file_name) as f:
        profiles_and_runs = [i.rstrip('\n').lower() for i in f.readlines()]
    fl_profiles_i = [n for n, l in enumerate(profiles_and_runs) if l.startswith('profile')]
    fl_runs_i = [n for n, l in enumerate(profiles_and_runs) if l.startswith('run')]
    fl_profiles = [(profiles_and_runs[n], daystr+profiles_and_runs[n+1], daystr+profiles_and_runs[n+2]) for n in fl_profiles_i]
    fl_runs = [(profiles_and_runs[n], daystr+profiles_and_runs[n+1], daystr+profiles_and_runs[n+2]) for n in fl_runs_i]
    res = fl_profiles+fl_runs
    if timesorted:
        from operator import itemgetter
        return sorted(res, key=itemgetter(2))
    else:
        return res
