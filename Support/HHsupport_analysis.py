#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 23:34:42 2022

@author: hhchoi1022
"""
#%%
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

from astropy.table import Table
from astropy.table import unique
from astropy.table import Column
from astropy.table import join
from astropy.table import vstack
try:
    import pyphot
    from pyphot import unit
except:
    pass


from scipy.interpolate import UnivariateSpline


#%%Visualization##################################################
##################################################################
##################################################################
##################################################################
##################################################################

def formatting_mosfit(mjd,
                       mag,
                       e_mag,
                       filter_,
                       observatory):
    
    formatted_tbl = Table()
    formatted_tbl['mjd'] = mjd
    formatted_tbl['mag'] = mag
    formatted_tbl['e_mag'] = e_mag
    formatted_tbl['filter'] = filter_
    formatted_tbl['telescope'] = observatory
    return formatted_tbl

def formatting_SNcosmo(mjd,
                       mag,
                       e_mag,
                       filter_,
                       magsys,
                       ):
    tbl = Table()
    tbl['mjd'] = mjd
    tbl['band'] = filter_
    tbl['flux'] = 10**((mag -25)/-2.5)
    tbl['fluxerr'] = e_mag*tbl['flux']*2.303/2.5
    tbl['zp'] = 25
    tbl['magsys'] = magsys
    tbl.sort('band')
    tbl.sort('mjd')
    return tbl

def load_marker_keys(length = None):
    '''
    parameters
    ----------
    {lenth}
    
    returns 
    -------
    1. marker_key with the length 
    
    notes 
    -----
    Order = ['o','s','^','P','*','D','p','+','X']
    -----
    '''
    marker = ['o','s','^','P','*','D','p','+','X']
    if length == None:
        return marker
    else:
        result_marker = marker[:length]
        return result_marker

def load_filt_keys(filter_key = None):
    '''
    parameters
    ----------
    the {list of filter}
    
    returns 
    -------
    1. color_key(matplotlib)
    2. offset_key(offset bw/ filters)
    3. filter_key(sncosmo)
    4. filter_key(pyphot)
    5. name_key(color_key+ offset_key)
    
    notes 
    -----
    5. name_key is for visualizing labels when plotting 
    -----
    '''
    color_key = dict(
                     U = 'cyan', 
                     B = 'b',
                     V = 'g',
                     R = 'r',
                     I = 'k',
                     u = 'darkcyan',
                     g = 'lightseagreen',
                     r = 'orange',
                     i = 'lightcoral',
                     z = 'gray'
                     )
    offset_key = dict(
                     U = -4, 
                     B = -2,
                     V = 0,
                     R = 1,
                     I = 3,
                     u = -5,
                     g = -1,
                     r = 2,
                     i = 4,
                     z = 5
                     )
    
    salt2_key = dict(
                     U = 'standard::u',
                     B = 'standard::b',
                     V = 'standard::v',
                     R = 'standard::r',
                     I = 'standard::i',
                     u = 'sdssu',
                     g = 'sdssg',
                     r = 'sdssr',
                     i = 'sdssi',
                     z = 'sdssz')
    
    pyphot_key = dict(
                     U = 'GROUND_JOHNSON_U',
                     B = 'GROUND_JOHNSON_B',
                     V = 'GROUND_JOHNSON_V',
                     R = 'GROUND_COUSINS_R',
                     I = 'GROUND_COUSINS_I',
                     u = 'SDSS_u',
                     g = 'SDSS_g',
                     r = 'SDSS_r',
                     i = 'SDSS_i',
                     z = 'SDSS_z')

    if filter_key ==None:    
        filter_key = [
                     'U',
                     'B',
                     'V',
                     'R',
                     'I',
                     'u',
                     'g',
                     'r',
                     'i',
                     'z'
                     ]
    result_color = {k: color_key[k] for k in filter_key if k in color_key}
    result_offset = {k: offset_key[k] for k in filter_key if k in offset_key}
    result_saltkey = {k: salt2_key[k] for k in filter_key if k in salt2_key}
    result_pyphotkey = {k: pyphot_key[k] for k in filter_key if k in pyphot_key}
    label_row = {k: k+'{:+}'.format(offset_key[k]) for k in filter_key if k in color_key}
    return result_color, result_offset, result_saltkey, result_pyphotkey, label_row

#%%Physical value operation#######################################
##################################################################
##################################################################
##################################################################
##################################################################

# Constants
c = 2.9979e10 # cm/s


def nu_to_lamb(nu):
    '''
    parameters
    ----------
    1. nu : frequency (Hz)
    
    returns 
    -------
    1. lamb : wavelength (AA)
    
    notes 
    -----
    lamb = c(speed of light) / nu
    -----
    '''
    lamb = c * 1e8 / nu
    return lamb

def lamb_to_nu(lamb):
    '''
    parameters
    ----------
    1. lamb : wavelength (AA)
    
    returns 
    -------
    1. nu : frequency (Hz)
    
    notes 
    -----
    nu = c(speed of light) / lamb
    -----
    '''
    nu = c * 1e8 / lamb
    return nu

def lflux_to_nuflux(lflux, wl_AA):
    '''
    parameters
    ----------
    1. f_lamb : flux [erg/s/cm^2/AA]
    2. Wavelength : [AA]
    
    returns 
    -------
    1. f_nu : flux [erg/s/cm^2/Hz]
    
    notes 
    -----
    -----
    '''
    f_nu = lflux  * wl_AA * wl_AA*1e-8 / c
    return f_nu

def nuflux_to_lflux(nuflux, wl_AA):
    '''
    parameters
    ----------
    1. nuflux [erg/s/cm^2/Hz]
    2. wl_AA : wavelength [AA]
    
    returns 
    -------
    1. f_lamb : flux [erg/s/cm^2/AA]
    
    notes 
    -----
    -----
    '''
    f_lamb = nuflux * c / wl_AA / wl_AA / 1e-8
    return f_lamb

def fnu_to_mag(f_nu):
    '''
    parameters
    ----------
    {f_nu[erg/s/cm^2/Hz]}
    
    returns 
    -------
    1. AB magnitude
    
    notes 
    -----
    -----
    '''
    mag = -2.5*np.log10(f_nu)-48.6
    return mag

def flamb_to_mag(f_lamb, wave_AA):

    '''
    parameters
    ----------
    {f_lamb[erg/s/cm^2/AA]}
    
    returns 
    -------
    1. AB magnitude
    
    notes 
    -----
    -----
    '''
    f_nu = nuflux_to_lflux(f_lamb, wave_AA)
    mag = -2.5*np.log10(f_nu)-48.6
    return mag

def mag_to_fnu(ABmag):
    '''
    parameters
    ----------
    1. AB magnitude
    
    returns 
    -------
    1. fnu [erg/s/cm^2/Hz]
    
    notes 
    -----
    -----
    '''
    fnu = 10**((ABmag + 48.6)/(-2.5)) 
    return fnu

def mag_to_flamb(ABmag,
                 wl):

    '''
    parameters
    ----------
    1. AB magnitude
    2. wavelength [AA]
    
    returns 
    -------
    1. flamb [erg/s/cm^2/Hz]
    
    notes 
    -----
    -----
    '''
    fnu = mag_to_fnu(ABmag)
    flamb = nuflux_to_lflux(fnu, wl)
    return flamb

def mag_to_flux(mag, zp = 25):
    '''
    parameters
    ----------
    {magnitude}
    
    returns 
    -------
    1. arbitrary flux with the zeropoint
    
    notes 
    -----
    -----

    '''
    flux = 10**((mag -zp)/-2.5)
    return flux

def flux_to_mag(flux, zp = 25):
    '''
    parameters
    ----------
    {flux} with the zeropoint
    
    returns 
    -------
    1. magnitude
    
    notes 
    -----
    -----
    '''
    mag = -2.5*np.log10(flux) + zp
    return mag

#%%Astropy Table Operation########################################
##################################################################
##################################################################
##################################################################
##################################################################

def remove_rows_table(tbl, column_key, remove_keys):
    '''
    parameters
    ----------
    1. {astropy table}
    2. {column name} in the table
    3. {remove_keys} to be removed from the values in the column of the table
    
    returns 
    -------
    1. table with the removed components 
    
    notes 
    -----
    This code removes rows matching {remove_keys}. 
    remove_keys can be either "string" or "list"
    -----
    '''
    if isinstance(remove_keys, str):
        remove_mask = tbl[column_key] == remove_keys
        remove_idx = np.where(remove_mask == True)
        tbl.remove_rows(remove_idx)
    else:
        for remove_key in remove_keys:
            remove_mask = tbl[column_key] == remove_key
            remove_idx = np.where(remove_mask == True)
            tbl.remove_rows(remove_idx)
    return tbl

def transpose_table(tab, id_col_name='ID'):
    '''
    parameters
    ----------
    {table} to be transposed with the {column names}
    
    returns 
    -------
    1. tranposed table
    
    notes 
    -----
    
    -----
    '''
    # contents of the first column of the old table provide column names for the new table
    # TBD: check for duplicates in new_colnames & resolve
    new_colnames=tuple(tab[tab.colnames[0]])
    # remaining columns of old table are row IDs for new table 
    new_rownames=tab.colnames[1:]
    # make a new, empty table
    tab_tansposed=Table(names=new_colnames)
    # add the columns of the old table as rows of the new table
    for r in new_rownames:
        tab_tansposed.add_row(tab[r])
    if id_col_name != '':
        # add the column headers of the old table as the id column of new table
        tab_tansposed.add_column(Column(new_rownames, name=id_col_name),index=0)
    return(tab_tansposed)

def read_Polin2019(filename):
    '''
    parameters
    ----------
    {filename} of Polin 2019 data to be read
    
    returns 
    -------
    1. data table

    notes 
    -----
    This is specialized to read the table from Polin 2019.
    -----
    '''
    header = ['phase', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']
    tbl = ascii.read(filename, format = 'fixed_width')
    tbl.rename_columns(tbl.colnames, header)
    return tbl

def read_Polin2019_spec(filename):
    '''
    parameters
    ----------
    {filename} of Polin 2019 spectrum data to read
    
    returns 
    -------
    1. data table
    2. phase

    notes 
    -----
    This is specialized to read the table from Polin 2019.
    -----
    '''

    import h5py
    from astropy.table import hstack
    f = h5py.File(filename, 'r')
    Lnu = f['Lnu'][:]    
    mu = f['mu'][0]
    nu = f['nu'][:]
    lamb = nu_to_lamb(nu)
    Llamb = nuflux_to_lflux(Lnu, lamb)
    time = f['time'][:]/86400
    Lnu = np.flip(Lnu, axis = 1)
    Llamb = np.flip(Llamb, axis = 1)
    nu = np.flip(nu)
    lamb = np.flip(lamb)
    Llamb /= 4*np.pi* (10*3.086e18)**2
    Lnu /= 4*np.pi* (10*3.086e18)**2
    return Lnu, Llamb, lamb, time, mu

def synphot_Polin2019_spec(filename,
                           filters = 'UBVRIugriz'):
    '''
    parameters
    ----------
    1. {filename} of Polin 2019 spectrum data to read
    2. filters : str = filterset for calculation
    
    returns 
    -------
    1. data table

    notes 
    -----
    This is specialized to read the table from Polin 2019.
    -----
    '''
    Lnu, Llamb, lamb, time, mu = read_Polin2019_spec(filename)
    spec = Llamb
    #spec = median_filter(spec, 10)
    filter_key = filters
    lib = pyphot.get_library()
    _, _, _, pyphot_key, _ = load_filt_keys(filter_key)
    model_phot = Table()
    model_phot['phase'] = time
    for filt_ in filter_key:
        filt_pyphot = lib[pyphot_key[filt_]]
        flux = filt_pyphot.get_flux(lamb*unit['AA'],spec*unit['ergs/s/cm**2/AA'], axis = 1)
        mag = -2.5*np.log10(flux.value) - filt_pyphot.AB_zero_mag
        model_phot[filt_] = mag
    return model_phot

def match_table(tbl1, tbl2, key, tolerance = 0.01):
    from astropy.table import hstack
    '''
    parameters
    ----------
    {two tables} to combine with the difference of the {key} smaller than the {tolerance}
    
    returns 
    -------
    1. combined table
    2. phase

    notes 
    -----
    Combined table have both columns of original tables. 
    They are horizontally combined in the order of tbl1, tbl2
    -----
    '''
    
    matched_tbl = Table()
    for obs in tbl1:
        ol_idx = (np.abs(obs[key] - tbl2[key]) < tolerance)
        if True in ol_idx:
            closest_idx = np.argmin(np.abs(obs[key]-tbl2[key]))
            compare_tbl = tbl2[closest_idx]
            compare_tbl = hstack([obs, compare_tbl])#join(obs, compare_tbl, keys = 'observatory', join_type = 'outer')
            matched_tbl = vstack([matched_tbl, compare_tbl])

    return matched_tbl

def binning_table(tbl, key, tolerance = 0.01):
    '''
    parameters
    ----------
    remove duplicates in the {table} with the difference of the {key} smaller than the {tolerance}
    
    returns 
    -------
    1. binned table 

    notes 
    -----
    
    -----
    '''

    rows = []
    for obs in tbl:
        ol_idx = (np.abs(obs[key] - tbl[key]) < tolerance)
        compare_tbl = tbl[ol_idx]
        row = []
        for val in compare_tbl.columns:
            if type(obs[val]) == np.str_:
                result_val = obs[val]
            elif type(obs[val]) ==np.ma.core.MaskedConstant:
                result_val = obs[val]
            else:
                result_val = round(np.mean(compare_tbl[val]),4)
            row.append(result_val)
        rows.append(row)
    binned_tbl = Table(data = np.array(rows),
                       dtype = compare_tbl.dtype,
                       names = compare_tbl.columns)
    result_tbl = unique(binned_tbl, key)
    return result_tbl

def groupping_table(tbl, key, tolerance = 0.1):
    '''
    parameters
    ----------
    group components in the {table} with the difference of the {key} smaller than the {tolerance}
    
    returns 
    -------
    1. groupped tables

    notes 
    -----
    
    -----
    '''
    i = 0
    table = tbl.copy()
    table['group'] = 0
    groupped_tbl = Table()
    while len(table) >= 1:
        group_idx = (np.abs(table[0][key] - table[key]) < tolerance)
        group_tbl = table[group_idx]
        group_tbl['group'] = i
        remove_idx = np.where(group_idx == True)
        table.remove_rows(remove_idx)
        groupped_tbl = vstack([group_tbl,groupped_tbl])
        i += 1
    
    return groupped_tbl



def read_HESMA(filename):
    '''
    parameters
    ----------
    HESMA {filename}
    
    returns 
    -------
    there are various type of files in HESMA database.
    
    [[[[Spectrum file]]]]
    1. wavelength
    2. phase[days]
    3. data table in flux_lambda [erg/s/cm^2/AA]
    4. data table in flux_nu [erg/s/cm^2/Hz]
    5. data table in AB magnitude
    6. synthetic photometry data table in AB magnitude 
    
    [[[[Bolometric light curve file]]]]
    1. phase[days]
    2. data table in flux [erg/s]

    [[[[early light curve file]]]]
    1. phase[days]
    2. data table in AB magnitude
    notes 
    -----
    
    -----
    '''      
    
                     
    def read_spectbl(filename,
                     filter_key ='UBVRIugriz'):
        tbl = ascii.read(filename)
        days = list(tbl[0])
        for i, colname in enumerate(tbl.copy().columns):
            tbl.rename_column(colname, days[i])
        tbl = tbl[1:]
        wl = list(tbl['0.0'])
        tbl.remove_column('0.0')
        f_lamb_tbl = tbl
        days = days[1:]
        f_nu_tbl = Table()
        mag_tbl = Table()
        for day in days:
            f_lamb = np.array(list(f_lamb_tbl[str(day)]))
            f_nu = lflux_to_nuflux(f_lamb, wl)
            mag = fnu_to_mag(f_nu)
            f_nu_tbl[f'{day}'] = f_nu
            mag_tbl[f'{day}'] = mag
        lib = pyphot.get_library()
        _, _, _, pyphot_key, _ = load_filt_keys(filter_key)
        synth_phot = Table()
        synth_phot['filter'] = list(filter_key)
        synth_phot.add_index('filter')
        for day in days:
            magset = []
            for filt_ in filter_key:
                filt_pyphot = lib[pyphot_key[filt_]]
                flux = filt_pyphot.get_flux(wl*unit['AA'],f_lamb_tbl[str(day)]*unit['ergs/s/cm**2/AA'], axis = 1)
                mag = -2.5*np.log10(flux.value) - filt_pyphot.AB_zero_mag
                magset.append(mag)
            synth_phot[f'{day}'] = magset
        return wl, days, f_lamb_tbl, f_nu_tbl, mag_tbl, synth_phot
    
    def read_lctbl(filename):
        result_tbl = ascii.read(filename)
        result_tbl.rename_column('col1','days')
        result_tbl.rename_column('col2','luminosity')
        days = result_tbl['days']
        return days, result_tbl
    
    def read_earlylctbl(filename):
        result_tbl = ascii.read(filename)
        return result_tbl
    
    if 'spectra' in filename:
        return read_spectbl(filename)
    elif 'lightcurve' in filename:
        if 'early' in filename:
            return read_earlylctbl(filename)
        else:
            return read_lctbl(filename)
    else:
        raise ValueError(f'{os.path.basename(filename)} cannot be interpreted')
#%%Data operation#################################################
##################################################################
##################################################################
##################################################################
##################################################################

def interpolate_spline(x, y, weight = None, k = 3, smooth = 0.05, show = False):
    '''
    parameters
    ----------
    Calculate interpolation with UnivariateSpline of the data{x}{y}
    weight : weight for each data point
    smooth : sum((w[i] * (y[i]-spl(x[i])))**2, axis=0) <= smooth
    show : show interpolated data 
    returns 
    -------
    1. {Spline function} that demonstrate data points 

    notes 
    -----
    -----
    '''

    s = UnivariateSpline(x,y, w = weight, s= smooth, k = k)
    fig = None
    if show:
        xgrid = np.arange(np.min(x),np.max(x), 0.005)
        fig = plt.figure(dpi =300)
        plt.gca().invert_yaxis()
        plt.scatter(x,y,marker= '+', s = 5, c = 'r', label  ='Raw data')
        plt.plot(xgrid, s(xgrid), c = 'k', linewidth = 1, label = 'Interpolation')
        plt.legend(loc = 0)
    return s, fig

def interpolate_linear(xdata, ydata, xrange_min, xrange_max, nx = 1000):
    '''
    parameters
    ----------
    1. xdata : np.array or list
            x data array for interpolation
    2. ydata : np.array or list
            y data array for interpolation
    3. xrange_min : float
            min limit of the interpolated output
    4. xrange_max : float
            max limit of the interpolated output
    5. nx : int
            the number of interpo;lated output        
            
    returns 
    -------
    1. output : list 
            the list of two arrays(interpolated_x, interpolated_y)

    notes 
    -----
    -----
    '''
    xgrid = np.linspace(xrange_min, xrange_max, nx)
    ygrid = np.interp(xgrid, xdata, ydata)
    return [xgrid, ygrid]
