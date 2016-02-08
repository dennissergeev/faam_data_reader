# -*- coding: utf-8 -*-
"""
Functions to read data from 2DS and CDP cloud probes
"""
import datetime
import netCDF4 as nc
import numpy as np

from . import utils


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


def read_cdp_nc(fname, time2datetime=True):
    """
    Read core FAAM **cloud** data from a NetCDF file.

    Args:
    -----
        fname: str, file name
    Kwargs:
    -------
        time2datetime: boolean, optional.
                       If True, convert time array to `datetime.datetime` objects.
                       Defaults to True.
    Returns:
    --------
        cdp_lwc_dens: liquid water content density [kg m :sup:`-3`]
        cdp_time: array-like of observations time

    """
    with nc.Dataset(fname) as cdp:
        cdp_time = cdp['Time']
        if time2datetime:
            if hasattr(cdp_time, 'units'):
                tbase, tstep_sec = utils.timestr2datetime(cdp_time.units)
                arr_sec2datetime = np.vectorize(lambda x: tbase + datetime.timedelta(seconds=int(x)*tstep_sec))
                cdp_time = arr_sec2datetime(cdp_time[:])
        else:
            cdp_time = cdp_time[:]

        ch_lims = np.vstack((cdp['CDP_D_L_NOM'][:], cdp['CDP_D_U_NOM'][:])) # Particle diameter lower and upper limits for each channel
        ch_mean_diam = np.mean(ch_lims,0) # Mean diameter for each channel
        ch_mean_vol = 4./3*np.pi*(0.5*ch_mean_diam)**3 # Mean volume for each channel

        h2o_d=999.97*1e3/(1e6)**3 # Water density in (g um-3)
        # Test (requires iris package):
        # a = iris.unit.Unit('kg m-3')
        # b = iris.unit.Unit('g um-3')
        # a.convert(999.97, b)
        # >>> 1e-12

        ch_mean_mass = ch_mean_vol*h2o_d # Mean mass for each channel
        assert len(ch_mean_mass) == len(cdp['CDP_D_L_NOM'][:]), 'Check the shape of arrays'

        cdp_lwc_dens_all_ch = []
        for ich, mass in enumerate(ch_mean_mass):
            cdp_conc = cdp['CDP_{0:02d}'.format(ich+1)][:] # droplet conc. in channel ich
            cdp_conc[cdp['CDP_FLAG'][:]!=0] = np.nan
            cdp_lwc_g_per_m3 = cdp_conc*mass*(1e2)**3
            cdp_lwc_dens_all_ch.append(cdp_lwc_g_per_m3*1e-3) # Append array of droplet densities in (kg m-3)
        cdp_lwc_dens = sum(np.array(cdp_lwc_dens_all_ch))

    return cdp_lwc_dens, cdp_time
