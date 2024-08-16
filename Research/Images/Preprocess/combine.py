#%%
from ccdproc import Combiner, CCDData
import numpy as np
from astropy.io import fits
from astropy.time import Time
import os

def bn_median(masked_array, 
              axis=None
              ):
    """

    parameters
    ----------
    masked_array : `numpy.ma.masked_array`
        Array of which to find the median.
    axis : optional, int 
        Axis along which to perform the median. Default is to find the median of
        the flattened array.
        
    returns
    ----------
    
    notes
    ----------
    Source code from Gregory S.H. Paek
    Perform fast median on masked array
    ----------
    """
    
    import numpy as np
    import bottleneck as bn
    data = masked_array.filled(fill_value=np.NaN)
    med = bn.nanmedian(data, axis=axis)
    # construct a masked array result, setting the mask from any NaN entries
    return np.ma.array(med, mask=np.isnan(med))

def combine_img(filelist,
                clip = 'sigma',
                combine = 'median',
                scale = 'zero',
                prefix = 'com_',
                refim = None,
                filename = None,
                
                # Clipping
                clip_sigma_low = 2,
                clip_sigma_high = 5,
                clip_minmax_min = 3,
                clip_minmax_max = 3,
                clip_extrema_nlow = 1,
                clip_extrema_nhigh = 1,
                
                # header 
                show_sum_exposure = True,
                show_ncombine = True,
                show_median_JD = True
                ):
    '''
    parameters
    ----------
    1. filelist : list or np.array 
            filelist to be combined
    2. clip : str
            method for clipping [None, minmax, sigma, extrema] (sigma)
    3. combine : str
            method for combining [mean, median, sum] (median)
    4. scale : str
            method for scaling [None, zero, multiply] (zero)
    5. prefix : str
            prefix of the combined image
    
    2.1. clip_sigma_low : optional, int
            Threshold for rejecting pixels that deviate below the baseline value.
    2.2. clip_sigma_high : optional, int
            Threshold for rejecting pixels that deviate above the baseline value.    
    2.3. clip_minmax_min : optional, int
            If not None, all pixels with values below min_clip will be masked.
    2.4. clip_minmax_max : optional, int
            If not None, all pixels with values above min_clip will be masked.
    2.5. clip_extrema_nlow : optional, int
            If not None, the number of low values to reject from the combination.
    2.6. clip_extrema_nhigh : optional, int
            If not None, the number of high values to reject from the combination.
    
    returns 
    -------
    1. outputname : str
            absolute path of the combined image
    
    notes 
    -----
    For more information : https://ccdproc.readthedocs.io/en/latest/image_combination.html
    -----
    '''
    print('Combining... \n')
    print(60*'=')    
    
    ccdlist = []
    hdr = fits.getheader(filelist[0])
    for file in filelist:
        ccd = CCDData.read(file, unit='adu')
        ccdlist.append(ccd)

    # Header
    if show_ncombine:
        hdr['NCOMBINE'] = int(len(filelist))
        for i, file in enumerate(filelist):
            hdr[f'COMBINE{i+1}'] = os.path.basename(file)
    if show_sum_exposure:
        hdr['EXPTIME'] = int(np.sum([inim.header['EXPTIME'] for inim in ccdlist]))
    if show_median_JD:
        if 'JD' in hdr.keys():
            hdr['JD'] = Time(np.mean([inim.header['JD'] for inim in ccdlist]), format='jd').value

    combiner = Combiner(ccdlist, dtype = np.float32)
    if scale == 'multiply':
        scaling_func = lambda arr: 1/np.ma.average(arr)
        combiner.scaling = scaling_func
    elif (scale == 'zero') & (refim == None):
        averages = [np.mean(ccddata) for ccddata in ccdlist]
        delvalues = averages - averages[0]
        for i, delvalue in enumerate(delvalues):
            ccdlist[i].data = ccdlist[i].data-delvalue
        combiner = Combiner(ccdlist, dtype = np.float32)

    # Clipping 
    if clip == 'minmax':
        combiner.minmax_clipping(min_clip=clip_minmax_min, max_clip=clip_minmax_max)
    if clip == 'sigma':
        combiner.sigma_clipping(low_thresh=clip_sigma_low, high_thresh=clip_sigma_high, func=np.ma.median)
    if clip == 'extrema':
        combiner.clip_extrema(nlow=clip_extrema_nlow, nhigh=clip_extrema_nhigh)
        
    # Combining 
    if combine == 'median':
        combined = combiner.median_combine(median_func = bn_median)
    if combine == 'mean':
        combined = combiner.average_combine()
    if combine == 'sum':
        combined = combiner.sum_combine()
    
    combined.header = hdr
    if filename == None:
        outputname = f'{os.path.dirname(filelist[0])}/{prefix}{os.path.basename(filelist[0])}'
    else:
        outputname = filename
    if (len(filelist) == 1):
            
        ccd.header = hdr
        ccd.write(outputname, overwrite=True, format = 'fits')
    else:
        combined.write(outputname, overwrite=True, format = 'fits')        
    print(f'Input image\n')
    print(60*'=')
    for file in filelist:
        print(file)
    print(60*'='+'\n')
    print(f'Combined path = {outputname}')
    
    return outputname


#%%
import glob
imlist = glob.glob('/home/hhchoi1022/Desktop/GRB221009A/KCT_STX16803/i/221013/align*.fit')

# %%
combine_img(imlist)
# %%
