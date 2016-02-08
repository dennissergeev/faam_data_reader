# -*- coding: utf-8 -*-
import datetime
import numpy as np

def uv2wdir(u, v):
    return 180.+180./np.pi*np.arctan2(u, v)


def uv2wspd(u, v):
    return np.sqrt(u**2 + v**2)


def timestr2datetime(tref, fmt='%Y-%m-%d%H:%M:%S'):
    """
    Parse a timestring and return a datetime object
    """
    spl = tref.rsplit()
    if spl[-1][0] == '+' or spl[-1] == 'UTC':
        date_and_time = spl[-3] + spl[-2]
    else:
        date_and_time = spl[-2] + spl[-1]
    dt = datetime.datetime.strptime(date_and_time, fmt)

    if   spl[0].lower() == 'days':
        nsec = 86400
    elif spl[0].lower() == 'hours':
        nsec = 3600
    elif spl[0].lower() == 'minutes':
        nsec = 60
    elif spl[0].lower() == 'seconds':
        nsec = 1
    else:
        nsec = 0
        print('Warning: unrecognized time unit, returned 0 for timestep')

    return dt, nsec


def filt_miss_row(*argin, **kwargs):
    """
    Filter N given numpy arrays by removing a row from all of the arrays
     if at least one array has a missing value in that row

    If a key word 'miss_value' is passed then it filters out the corresponding values.
    Otherwise, if input arrays have 'fill_value' property, fill_value values are removed,
    Otherwise, NaNs are filtered out.
    """
    nargs = len(argin)

    if nargs > 0:

        args = [x.flat[:] for x in argin]
        n = args[0].size

        if not all(x.size == n for x in args):
            raise NameError('Vectors must the same length')

        work = np.array(argin).T
        if 'miss_val' in kwargs:
            return list(work[np.all(work != kwargs['miss_val'], axis=1)].T) 
        else:
            if hasattr(argin[0],'fill_value'):
                return list(work[np.all(work != argin[0].fill_value, axis=1)].T) 
            else:
                return list(work[np.all(~np.isnan(work), axis=1)].T) 
