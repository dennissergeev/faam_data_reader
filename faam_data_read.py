# -*- coding: utf-8 -*-
"""
Functions to read FAAM data:
 * core processed data
 * dropsonde data
 * 2DS and CDP instruments data
"""
import datetime
import netCDF4 as nc
import numpy as np

import var_utils as var
import phys_meteo as met

class ObsData:
    """Generic class for storing several fields of observational data.

    Contains methods `__init__` and `__call__` for initialising and adding
    fields to this class.
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
    def __call__(self,**kwds):
        self.__dict__.update(kwds)


class DsFld:
    """A class for storing dropsonde data.

    Contains several attributes:
        raw: array-like, 'raw' data
        fil: array-like, 'filtered' data
        units: string, units
        long_name: string, name
    """
    def __init__(self, raw=np.array([]), fil=np.array([]), units='', long_name=''):
        self.raw = raw
        self.fil = fil
        self.units = units
        self.long_name = long_name


class FaamFld:
    """A class for storing core FAAM aircraft data.

    Contains several attributes:
        val: array-like, data values
        units: str, units
        long_name: str, name
    """
    def __init__(self, val=np.array([]), units='', long_name=''):
        self.val = val
        self.units = units
        self.long_name = long_name


def read_ds_nc(fname, flds=None, time2datetime=False):
    """Read dropsonde data from a NetCDF file.

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
    fil_list = var.filt_miss_row(*[getattr(dum,ii).raw for ii in flds_list])

    data = ObsData()
    for i, j in enumerate(fil_list):
        data(**{flds_list[i]:DsFld(raw=getattr(dum,flds_list[i]).raw,\
                                   fil=j,\
                                   units=getattr(dum,flds_list[i]).units,\
                                   long_name=getattr(dum,flds_list[i]).long_name)})

    if time2datetime and 'time' in flds:
        if hasattr(data.time, 'units'):
            tbase, tstep_sec = var.timestr2datetime(data.time.units)
            arr_sec2datetime = np.vectorize(lambda x: tbase + datetime.timedelta(seconds=x*tstep_sec))
            data.time.fil = arr_sec2datetime(data.time.fil)

    return data


def read_faam_nc(fname, flds=None, time2datetime=False, calc_wspd=True, calc_wdir=True):
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

        data(**{i:FaamFld(val=ncdata[:].squeeze(),units=ncfld.units,long_name=ncfld.long_name)})

    if hasattr(data,'u') and hasattr(data,'v'):
        if calc_wspd:
            data.wspd = FaamFld(var.uv2wspd(data.u.val,data.v.val),data.u.units,'wind speed derived from aircraft instruments and GIN')
        if calc_wdir:
            data.wdir = FaamFld(var.uv2wdir(data.u.val,data.v.val),'deg','wind direction')

    if time2datetime and 'time' in flds:
        if hasattr(data.time, 'units'):
            tbase, tstep_sec = var.timestr2datetime(data.time.units)
            arr_sec2datetime = np.vectorize(lambda x: tbase + datetime.timedelta(seconds=int(x)*tstep_sec))
            data.time.val = arr_sec2datetime(data.time.val)

    return data


def read_2ds_hdf(fname, tbase=datetime.datetime(2013, 3, 26, 0, 0, 0), time2datetime=True):
    """Read 2DS data from a HDF5 file. Requires `h5py` package.

    Sum data from all channels and convert densities to [kg m :sup:`-3`].

    Args:
    -----
        fname: str, file name
    Kwargs:
    -------
        tbase: datetime.datetime, time start. Defaults to datetime.datetime(2013, 3, 26, 0, 0, 0).
        time2datetime: boolean, optional.
                       If True, convert time array to `datetime.datetime` objects
                       adding time values to `tbase` kwarg.
                       Defaults to True.
    Returns:
    --------
        time: array-like of observations time
        fwc: frozen water content density [kg m :sup:`-3`]
        lwc: liquid water content density [kg m :sup:`-3`]
        other: small particles density [kg m :sup:`-3`]

    """

    import h5py

    ds2 = h5py.File(fname)
    time = [ds2[i] for i in ds2.keys() if 'time' in i.lower()][0]
    if time2datetime:
        time = np.array([tbase + datetime.timedelta(seconds=i) for i in time.value])
    else:
        time = time.value

    fwc = np.nansum([ds2[i] for i in ds2.keys() if 'I_MD' in i][0].value,1) # Frozen water content density
    lwc = np.nansum([ds2[i] for i in ds2.keys() if 'R_MD' in i][0].value,1) # Liquid water content density
    other = fwc = np.nansum([ds2[i] for i in ds2.keys() if 'S_MD' in i][0].value,1)  # small particles density

    fwc, lwc, other = [i*1e-3 for i in (fwc, lwc, other)] # Convert to kg m-3

    return time, fwc, lwc, other

def read_faam_cdp_nc(fname, time2datetime=True, use='xray'):
    """Read core FAAM **cloud** data from a NetCDF file. Requires `xray` package.

    Args:
    -----
        fname: str, file name
    Kwargs:
    -------
        time2datetime: boolean, optional.
                       If True, convert time array to `datetime.datetime` objects.
                       Defaults to True.
        use: str, optional
             Defaults to 'xray'.
    Returns:
    --------
        cdp_lwc_dens: liquid water content density [kg m :sup:`-3`]
        cdp_time: array-like of observations time

    TODO: netCDF interface
    """
    if use.lower() == 'xray':
        import xray
        cdp = xray.open_dataset(fname)
        cdp_time = cdp.Time.values
        if time2datetime:
            cdp_time = cdp_time.astype('<M8[us]').astype(datetime.datetime) # Time as datetime objects

        ch_lims = np.vstack((cdp.CDP_D_L_NOM.values, cdp.CDP_D_U_NOM.values)) # Particle diameter lower and upper limits for each channel
        ch_mean_diam = np.mean(ch_lims,1) # Mean diameter for each channel
        ch_mean_vol = 4./3*np.pi*(0.5*ch_mean_diam)**3 # Mean volume for each channel

        h2o_d=999.97*1e3/(1e6)**3 # Water density in (g um-3)
        # Test (requires iris package):
        # a = iris.unit.Unit('kg m-3')
        # b = iris.unit.Unit('g um-3')
        # a.convert(999.97, b)
        # >>> 1e-12

        ch_mean_mass = ch_mean_vol*h2o_d # Mean mass for each channel

        cdp_lwc_dens_all_ch = []
        for ich, mass in enumerate(ch_mean_mass):
            cdp_conc = cdp['CDP_{0:02d}'.format(ich+1)].values # droplet conc. in channel ich
            cdp_conc[cdp.CDP_FLAG.values!=0] = np.nan
            cdp_lwc_dens_all_ch.append(cdp_conc*mass*(1e2)**3/1e3) # Append array of droplet densities in (kg m-3)
        cdp_lwc_dens = sum(np.array(cdp_lwc_dens_all_ch))
    else:
        raise NotImplementedError('Only xray interface works now')

    return cdp_lwc_dens, cdp_time


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
    profiles_and_runs = [i.rstrip('\n').lower() for i in open(text_file_name).readlines()]
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
