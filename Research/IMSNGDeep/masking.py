#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 15:58:05 2022

@author: hhchoi1022
"""
#%%
import astropy.io.fits as fits
import numpy as np
from astropy.stats import SigmaClip
from photutils.background import SExtractorBackground
import os

image = '/data2/IMSNGDEEP/Analysis/masking/NGC3147.com.fits'

detthresh = 20

conf_param = dict(
    # Default configuration file for Source Extractor 2.25.0
    # EB 2018-02-08
    #
     
    #-------------------------------- Catalog ------------------------------------
     
    CATALOG_NAME     ='mask.cat',       # name of the output catalog
    CATALOG_TYPE     ='ASCII_HEAD',     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
                                    # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC
    PARAMETERS_NAME  ='mask.SEparam',  # name of the file containing catalog contents
     
    #------------------------------- Extraction ----------------------------------
     
    DETECT_TYPE      ='CCD',            # CCD (linear) or PHOTO (with gamma correction)
    DETECT_MINAREA   =5,              # min. # of pixels above threshold
    DETECT_MAXAREA   =0,              # max. # of pixels above threshold (0=unlimited)
    THRESH_TYPE      ='RELATIVE',       # threshold type: RELATIVE (in sigmas)
                                    # or ABSOLUTE (in ADUs)
    DETECT_THRESH    =1.5,            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
    ANALYSIS_THRESH  =1.5,            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
     
    FILTER           ='Y',              # apply filter for detection (Y or N)?
    FILTER_NAME      ='default.conv',   # name of the file containing the filter
     
    DEBLEND_NTHRESH  =32,             # Number of deblending sub-thresholds
    DEBLEND_MINCONT  =0.005,          # Minimum contrast parameter for deblending
     
    CLEAN            ='Y',              # Clean spurious detections? (Y or N)?
    CLEAN_PARAM      =1.0,            # Cleaning efficiency
     
    MASK_TYPE        ='CORRECT',        # type of detection MASKing: can be one of
                                    # NONE, BLANK or CORRECT
     
    #-------------------------------- WEIGHTing ----------------------------------
    
    WEIGHT_TYPE      ='NONE',           # type of WEIGHTing: NONE, BACKGROUND,
                                    # MAP_RMS, MAP_VAR or MAP_WEIGHT
    RESCALE_WEIGHTS  ='Y',              # Rescale input weights/variances (Y/N)?
    WEIGHT_IMAGE     ='weight.fits',    # weight-map filename
    WEIGHT_GAIN      ='Y',              # modulate gain (E/ADU) with weights? (Y/N)
    
    #------------------------------ Photometry -----------------------------------
     
    PHOT_APERTURES   ='5',              # MAG_APER aperture diameter(s) in pixels
    PHOT_AUTOPARAMS  ='2.5,3.5',       # MAG_AUTO parameters: <Kron_fact>,<min_radius>
    PHOT_PETROPARAMS ='2.0,3.5',       # MAG_PETRO parameters: <Petrosian_fact>,
                                    # <min_radius>
    PHOT_AUTOAPERS   ='0.0,0.0',        # <estimation>,<measurement> minimum apertures
                                    # for MAG_AUTO and MAG_PETRO
    PHOT_FLUXFRAC    =0.5,            # flux fraction[s] used for FLUX_RADIUS
     
    SATUR_LEVEL      =50000.0,        # level (in ADUs) at which arises saturation
    SATUR_KEY        ='SATURATE',       # keyword for saturation level (in ADUs)
     
    MAG_ZEROPOINT    =0.0,            # magnitude zero-point
    MAG_GAMMA        =4.0,            # gamma of emulsion (for photographic scans)
    GAIN             =0.0,            # detector gain in e-/ADU
    GAIN_KEY         ='GAIN',           # keyword for detector gain in e-/ADU
    PIXEL_SCALE      =1.0,            # size of pixel in arcsec (0=use FITS WCS info)
     
    #------------------------- Star/Galaxy Separation ----------------------------
     
    SEEING_FWHM      =3.5,            # stellar FWHM in arcsec
    STARNNW_NAME     ='default.nnw',    # Neural-Network_Weight table filename
     
    #------------------------------ Background -----------------------------------
     
    BACK_TYPE        ='AUTO',           # AUTO or MANUAL
    BACK_VALUE       =0.0,            # Default background value in MANUAL mode
    BACK_SIZE        =64,             # Background mesh: <size> or <width>,<height>
    BACK_FILTERSIZE  =3,              # Background filter: <size> or <width>,<height>
     
    BACKPHOTO_TYPE   ='LOCAL',         # can be GLOBAL or LOCAL
    BACKPHOTO_THICK  =24,             # thickness of the background LOCAL annulus
    BACK_FILTTHRESH  =0.0,            # Threshold above which the background-
                                    # map filter operates
     
    #------------------------------ Check Image ----------------------------------
     
    CHECKIMAGE_TYPE  ='SEGMENTATION',           # can be NONE, BACKGROUND, BACKGROUND_RMS,
                                    # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                    # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                    # or APERTURES
    CHECKIMAGE_NAME  ='mask.fits',     # Filename for the check-image
     
    #--------------------- Memory (change with caution!) -------------------------
     
    MEMORY_OBJSTACK  =3000,           # number of objects in stack
    MEMORY_PIXSTACK  =300000,         # number of pixels in stack
    MEMORY_BUFSIZE   =1024,           # number of lines in buffer
     
    #--------------------------- Experimental Stuff -----------------------------
    
    PSF_NAME         ='default.psf',    # File containing the PSF model
    PSF_NMAX         =1,              # Max.number of PSFs fitted simultaneously
    PATTERN_TYPE     ='RINGS-HARMONIC', # can RINGS-QUADPOLE, RINGS-OCTOPOLE,
                                    # RINGS-HARMONICS or GAUSS-LAGUERRE
    SOM_NAME         ='default.som'    # File containing Self-Organizing Map weights
    )
from astropy.wcs import WCS
from astropy.io import ascii
hdu = fits.open(image)[0]
hdr = hdu.header
wcs = WCS(hdr)

conf_param['DETECT_THRESH'] = detthresh


config = ''
for param in conf_param.keys():
    config += f'-{param} {conf_param[param]} '

os.chdir('/data2/sextractor')
os.system(f'source-extractor {image} {config}   ')

setable = ascii.read('/data2/sextractor/mask.cat')
setable['SKYSIG'] = setable['THRESHOLD']/detthresh
setable['SNR'] = 1/setable['MAGERR_AUTO']
#%%
# Objects cutouts
def CUTOUT(image, x_image, y_image, peeing, bkgvalue, factor_size = 15, cutoutsize = None, write=False, imname = None):
    from astropy.nddata import Cutout2D
    import os
    from astropy.wcs import WCS
    from astropy.io import fits
    if type(image) == str:
        hdu = fits.open(image)[0]
    else :
        hdu = image
    wcs = WCS(hdu.header)
    size = float(peeing * factor_size)
    if cutoutsize != None:
        size = cutoutsize
    position = [x_image, y_image]
    cutout = Cutout2D(hdu.data, position , size, wcs= wcs)
    hdu.data = cutout.data-bkgvalue
    hdu.header.update(cutout.wcs.to_header())
    curpath = os.path.dirname(image)
    if write == True:
        if imname == None:
            os.makedirs(f'{curpath}/cutout', exist_ok = True)
            os.chdir(f'{curpath}/cutout')
            hdu.writeto(f'{curpath}/cutout/{os.path.basename(image)}',overwrite=True)
        else:
            os.makedirs(f'{curpath}/cutout', exist_ok = True)
            os.chdir(f'{curpath}/cutout')
            hdu.writeto(f'{imname}', overwrite = True)
        return hdu
    else:
        return hdu

condition = [(setable['CLASS_STAR']>0.5)
             & (setable['FLAGS'] != 4)
             & (setable['SNR'] > 5000)]

setable = setable[condition]
x_image = setable[0]['X_IMAGE']
y_image = setable[0]['Y_IMAGE']
peeing = setable[0]['FWHM_IMAGE']

for star in setable:
    x_image = star['X_IMAGE']
    y_image = star['Y_IMAGE']
    peeing = star['FWHM_IMAGE']
    bkgvalue = star['BACKGROUND']
    number = star['NUMBER']
    CUTOUT(image, x_image, y_image, peeing, bkgvalue, write = True, imname = number)

#%%
# Gaussian fitting
setable


import glob
from scipy.optimize import curve_fit
for cutout_image in glob.glob('/data2/IMSNGDEEP/Analysis/masking/cutout/*'):
    
    index = int(os.path.basename(cutout_image))
    index = 85
    cutout_image = f'/data2/IMSNGDEEP/Analysis/masking/cutout/{index}'
    star = setable[setable['NUMBER'] == index]
    # Getting the 1D data
    data = fits.getdata(cutout_image)
    data1D = data[int(data.shape[0]/2)-1]
    
    def func(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))
    
    amplitude = np.max(data1D)
    sigma = star['FWHM_IMAGE']/2.35
    x = np.linspace(0, data.shape[0], 100)
    y = func(x, amplitude, int(data.shape[0]/2)-1, sigma)
    import matplotlib.pyplot as plt
    plt.title(index)
    plt.xlim(10,20)
    plt.ylim(0,1)
    plt.plot(data1D/star['THRESHOLD']*20)
    plt.plot(x,y/star['THRESHOLD']*20, c = 'k')
    plt.show()
    
    
#%%

bkgimage = '/data2/sextractor/mask.fits'
data = fits.getdata(bkgimage)
np.median(data)
import astropy.io.fits as fits
from astropy.stats import sigma_clipped_stats
import numpy as np
from astropy.stats import SigmaClip
from photutils.background import SExtractorBackground
import os
import parmap
from photutils.segmentation import make_source_mask
sigma_clip = SigmaClip(sigma=3.0)
bkg = SExtractorBackground(sigma_clip)

def mask_bg_norm(im):
	print(im)
	data=fits.getdata(im)
	fn=os.path.splitext(im)[0]
	bkg_value = bkg.calc_background(data)
	normdata=data/bkg_value
	mask = make_source_mask(data, nsigma=1, npixels=2)
	#mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
	#mean, median, std = sigma_clipped_stats(data, sigma=3.0)
	masked_data = ~mask * normdata
	fits.writeto( fn+'_mask.fits', masked_data,overwrite=True)

#%%
data=fits.getdata(bkgimage)
fn=os.path.splitext(bkgimage)[0]
bkg_value = bkg.calc_background(data)
normdata=data/bkg_value
mask = make_source_mask(data, nsigma=3, npixels = 5 )
masked_data = ~mask * normdata
fits.writeto( fn+'_mask.fits', masked_data,overwrite=True)
img_data = fits.getdata(image)
masked_imgdata = ~mask * img_data
fn=os.path.splitext(image)[0]
fits.writeto( fn+'_masked.fits', masked_imgdata,overwrite=True)
masked_image = '/data2/IMSNGDEEP/Analysis/masking/NGC3147.com_masked.fits'
mask_bg_norm(masked_image)
coverage_mask = (masked_data == 0)
coverage_mask



#%%






from astropy.stats import SigmaClip
from photutils.background import SExtractorBackground
import os, glob

image = '/data2/IMSNGDEEP/Analysis/masking/NGC1566.com.fits'

sexpath  ='/data2/sextractor'

conf_param = dict(
    # Default configuration file for Source Extractor 2.25.0
    # EB 2018-02-08
    #
    #-------------------------------- Catalog ------------------------------------
     
    CATALOG_NAME     ='mask.cat',       # name of the output catalog
    CATALOG_TYPE     ='ASCII_HEAD',     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
                                    # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC
    PARAMETERS_NAME  ='mask.SEparam',  # name of the file containing catalog contents
     
    #------------------------------- Extraction ----------------------------------
     
    DETECT_TYPE      ='CCD',            # CCD (linear) or PHOTO (with gamma correction)
    DETECT_MINAREA   =5,              # min. # of pixels above threshold
    DETECT_MAXAREA   =0,              # max. # of pixels above threshold (0=unlimited)
    THRESH_TYPE      ='RELATIVE',       # threshold type: RELATIVE (in sigmas)
                                    # or ABSOLUTE (in ADUs)
    DETECT_THRESH    =1.5,            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
    ANALYSIS_THRESH  =1.5,            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
     
    FILTER           ='Y',              # apply filter for detection (Y or N)?
    FILTER_NAME      ='default.conv',   # name of the file containing the filter
     
    DEBLEND_NTHRESH  =32,             # Number of deblending sub-thresholds
    DEBLEND_MINCONT  =0.005,          # Minimum contrast parameter for deblending
     
    CLEAN            ='Y',              # Clean spurious detections? (Y or N)?
    CLEAN_PARAM      =1.0,            # Cleaning efficiency
     
    MASK_TYPE        ='CORRECT',        # type of detection MASKing: can be one of
                                    # NONE, BLANK or CORRECT
     
    #-------------------------------- WEIGHTing ----------------------------------
    
    WEIGHT_TYPE      ='NONE',           # type of WEIGHTing: NONE, BACKGROUND,
                                    # MAP_RMS, MAP_VAR or MAP_WEIGHT
    RESCALE_WEIGHTS  ='Y',              # Rescale input weights/variances (Y/N)?
    WEIGHT_IMAGE     ='weight.fits',    # weight-map filename
    WEIGHT_GAIN      ='Y',              # modulate gain (E/ADU) with weights? (Y/N)
    
    #------------------------------ Photometry -----------------------------------
     
    PHOT_APERTURES   ='5',              # MAG_APER aperture diameter(s) in pixels
    PHOT_AUTOPARAMS  ='2.5,3.5',       # MAG_AUTO parameters: <Kron_fact>,<min_radius>
    PHOT_PETROPARAMS ='2.0,3.5',       # MAG_PETRO parameters: <Petrosian_fact>,
                                    # <min_radius>
    PHOT_AUTOAPERS   ='0.0,0.0',        # <estimation>,<measurement> minimum apertures
                                    # for MAG_AUTO and MAG_PETRO
    PHOT_FLUXFRAC    =0.5,            # flux fraction[s] used for FLUX_RADIUS
     
    SATUR_LEVEL      =50000.0,        # level (in ADUs) at which arises saturation
    SATUR_KEY        ='SATURATE',       # keyword for saturation level (in ADUs)
     
    MAG_ZEROPOINT    =0.0,            # magnitude zero-point
    MAG_GAMMA        =4.0,            # gamma of emulsion (for photographic scans)
    GAIN             =0.0,            # detector gain in e-/ADU
    GAIN_KEY         ='GAIN',           # keyword for detector gain in e-/ADU
    PIXEL_SCALE      =1.0,            # size of pixel in arcsec (0=use FITS WCS info)
     
    #------------------------- Star/Galaxy Separation ----------------------------
     
    SEEING_FWHM      =3.5,            # stellar FWHM in arcsec
    STARNNW_NAME     ='default.nnw',    # Neural-Network_Weight table filename
     
    #------------------------------ Background -----------------------------------
     
    BACK_TYPE        ='AUTO',           # AUTO or MANUAL
    BACK_VALUE       =0.0,            # Default background value in MANUAL mode
    BACK_SIZE        =64,             # Background mesh: <size> or <width>,<height>
    BACK_FILTERSIZE  =3,              # Background filter: <size> or <width>,<height>
     
    BACKPHOTO_TYPE   ='LOCAL',         # can be GLOBAL or LOCAL
    BACKPHOTO_THICK  =24,             # thickness of the background LOCAL annulus
    BACK_FILTTHRESH  =0.0,            # Threshold above which the background-
                                    # map filter operates
     
    #------------------------------ Check Image ----------------------------------
     
    CHECKIMAGE_TYPE  ='NONE',           # can be NONE, BACKGROUND, BACKGROUND_RMS,
                                    # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                    # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                    # or APERTURES
    CHECKIMAGE_NAME  ='mask.fits',     # Filename for the check-image
     
    #--------------------- Memory (change with caution!) -------------------------
     
    MEMORY_OBJSTACK  =3000,           # number of objects in stack
    MEMORY_PIXSTACK  =300000,         # number of pixels in stack
    MEMORY_BUFSIZE   =1024,           # number of lines in buffer
     
    #--------------------------- Experimental Stuff -----------------------------
    
    PSF_NAME         ='default.psf',    # File containing the PSF model
    PSF_NMAX         =1,              # Max.number of PSFs fitted simultaneously
    PATTERN_TYPE     ='RINGS-HARMONIC', # can RINGS-QUADPOLE, RINGS-OCTOPOLE,
                                    # RINGS-HARMONICS or GAUSS-LAGUERRE
    SOM_NAME         ='default.som'    # File containing Self-Organizing Map weights
    )
from astropy.wcs import WCS
from astropy.io import ascii

conf_param['CHECKIMAGE_TYPE'] = 'BACKGROUND'
conf_param['CHECKIMAGE_NAME'] = 'mask_bkg.fits'


config = ''
for param in conf_param.keys():
    config += f'-{param} {conf_param[param]} '

os.chdir(sexpath)
os.system(f'source-extractor {image} {config}   ')

bkgpath = f'{sexpath}/{conf_param["CHECKIMAGE_NAME"]}'

# 1st masking with SExtractor background map
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.segmentation import make_source_mask
def mask_bg_norm(image, sigfactor = 3):
    print(image)
    sigma_clip = SigmaClip(sigma=3.0)
    bkg = SExtractorBackground(sigma_clip)
    data=fits.getdata(image)
    fn=os.path.splitext(image)[0]
    bkg_value = bkg.calc_background(data)
    normdata=data/bkg_value
    #mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    mask = make_source_mask(data, nsigma=std*sigfactor, npixels=5)
    return ~mask, bkg_value

mask1, bkg_value = mask_bg_norm(bkgpath, sigfactor = 5)
mask1_int = mask1.astype(int)
img_data = fits.getdata(image)
img_hdr = fits.getheader(image)
masked1_img_data = img_data * mask1
fn=os.path.splitext(image)[0]
mask1_path  = fn + '.mask1.fits'
masked1_path  = fn + '.masked1.fits'
fits.writeto( mask1_path, mask1_int, img_hdr, overwrite=True)
fits.writeto( masked1_path, masked1_img_data, img_hdr, overwrite=True)

#
conf_param['BACK_TYPE'] = 'MANUAL'
conf_param['BACK_VALUE'] = round(float(bkg_value),2)
conf_param['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
conf_param['CHECKIMAGE_NAME'] = 'mask_tmp.fits'
conf_param['DETECT_THRESH'] = 1.5

config = ''
for param in conf_param.keys():
    config += f'-{param} {conf_param[param]} '

os.chdir(sexpath)
os.system(f'source-extractor {masked1_path} {config}   ')
mask2_data = fits.getdata(f'{sexpath}/{conf_param["CHECKIMAGE_NAME"]}')
mask2 = (mask2_data != 0)
mask2_int = mask2.astype(int)
img_data = fits.getdata(image)
img_hdr = fits.getheader(image)
masked2_img_data = img_data * mask2
fn=os.path.splitext(image)[0]
mask2_path  = fn + '.mask2.fits'
masked2_path  = fn + '.masked2.fits'
fits.writeto( mask2_path, mask2_int, img_hdr, overwrite=True)
fits.writeto( masked2_path, masked2_img_data, img_hdr, overwrite=True)

master_mask = ~(~mask1 | mask2)
master_mask = master_mask.astype(int)
masked_img_data = img_data * master_mask

fits.writeto(f'{fn}.mask.fits', master_mask, img_hdr, overwrite  =True)
fits.writeto(f'{fn}.masked.fits', masked_img_data, img_hdr, overwrite  =True)




