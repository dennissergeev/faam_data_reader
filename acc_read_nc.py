# -*- coding: utf-8 -*-
"""
@author: D. E. Sergeev
"""
import datetime
import netCDF4 as nc
import numpy as np
# My modules
import var_utils as var
from phys_meteo import uv2wspd, uv2wdir

class ObsData:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
    def __call__(self,**kwds):
        self.__dict__.update(kwds)

class DsFld:
    def __init__(self, raw=np.array([]), fil=np.array([]), units='', long_name=''):
        self.raw = raw
        self.fil = fil
        self.units = units
        self.long_name = long_name

class FaamFld:
    def __init__(self, val=np.array([]), units='', long_name=''):
        self.val = val
        self.units = units
        self.long_name = long_name

def read_ds_nc(fname, flds=None, time2datetime=False):
    """Reads several variables from a dropsonde processed netcdf file
    Stores variables in structure array 'data', instead of using NetCDF4 variable class
    Each variable in the 'data' array usually has the following fields
     - 'raw' (numbers themselves)
     - 'units'
     - 'long_name'
     - 'fil' values without any missing values in all columns
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
    fil_list = var.filt_miss_row(*[getattr(dum,ii).raw for ii in flds_list])

    data = ObsData()
    for i,j in enumerate(fil_list):
        data(**{flds_list[i]:DsFld(raw=getattr(dum,flds_list[i]).raw,\
                                   fil=j,\
                                   units=getattr(dum,flds_list[i]).units,\
                                   long_name=getattr(dum,flds_list[i]).long_name)})
    
    if time2datetime and 'time' in flds:
        if hasattr(data.time, 'units'):
            tbase, tstep_sec = var.timestr2datetime(data.time.units)
            arr_sec2datetime = np.vectorize(lambda x: tbase + datetime.timedelta(seconds=x*tstep_sec))                       
            dt = arr_sec2datetime(data.time.fil)
            data.time.fil = dt
    
    return data

def read_faam_nc(fname, flds=None, time2datetime=False):
    """Reads core FAAM data from a netcdf file
    Stores variables in structure array 'data'
    Each variable in the 'data' array usually has the following fields
     - 'raw' (numbers themselves)
     - 'units' (as in .nc file)
     - 'long_name' (as in .nc file)
    """

    if flds == None:
        flds = dict(time='Time',hgt='HGT_RADR',lon='LON_GPS',lat='LAT_GPS',ang='TRCK_GIN',\
                    uturb='U_C',vturb='V_C',wturb='W_C',\
                    u='U_NOTURB',v='V_NOTURB',\
                    pres='PS_RVSM',tdry='TAT_ND_R')

    f = nc.Dataset(fname)

    data = ObsData()
    for i in flds:
        ncfld = f.variables[flds[i]]
        ncdata = ncfld[:]
        try:
            ncflag = f.variables[flds[i]+'_FLAG']
            ncdata[ncflag[:]!=0] = float('nan')
        except KeyError:
            pass

        data(**{i:FaamFld(val=ncdata[:],units=ncfld.units,long_name=ncfld.long_name)})

    if hasattr(data,'u') and hasattr(data,'v'):
        data.wspd = FaamFld(uv2wspd(data.u.val,data.v.val),data.u.units,'wind speed derived from aircraft instruments and GIN')
        data.wdir = FaamFld(uv2wdir(data.u.val,data.v.val),'deg','wind direction')
        
    if time2datetime and 'time' in flds:
        if hasattr(data.time, 'units'):
            tbase, tstep_sec = var.timestr2datetime(data.time.units)
            arr_sec2datetime = np.vectorize(lambda x: tbase + datetime.timedelta(seconds=x*tstep_sec))
            dt = arr_sec2datetime(data.time.val)
            data.time.val = dt
            
    return data
