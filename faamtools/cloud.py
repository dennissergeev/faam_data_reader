# -*- coding: utf-8 -*-
"""
Functions to read data from 2DS and CDP cloud probes
"""
import datetime
import numpy as np

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
