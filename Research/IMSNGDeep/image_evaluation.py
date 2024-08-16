#!/usr/bin/env python
# coding: utf-8

# This code is written by Hyeonho Choi at 2022/02/11
# Some parts are refered from the photometry code by Gregory Paek gppy2
# To use this code, conversion.py is needed for magnitude conversion between different filterset.
# Depth, Seeing, ZP... will be updated after running this code (UPDATE IMAGE)

#%%
import os
os.chdir('/home/hhchoi1022/Desktop/Gitrepo/Reference_make')
from Observation import cross_match
from conversion import PANSTARRS1_to_SDSS
from conversion import PANSTARRS1_to_JH





def TARGETSIZE_PIXEL(target, obs_info, observatory):
    from astropy.io import ascii
    
    IMSNG_fieldlist =  ascii.read(f'/home/hhchoi1022/Desktop/Gitrepo/config/alltarget.dat')
    pixsize = obs_info['pixelscale']
    
    targetsize = IMSNG_fieldlist[IMSNG_fieldlist['obj'] == target]['maxaxis']*60/pixsize
    if targetsize <=64:
        bkgsize = 64
    elif targetsize > 64 and targetsize <=128:
        bkgsize = 128
    elif targetsize >128 and targetsize <=400:
        bkgsize = 256  
    elif targetsize >400:
        bkgsize = 512  
    return bkgsize

def LOAD_APASS(target, filter_):
    from conversion import SDSS_to_JH
    from astropy.io import ascii
    from astropy.table import Table
    
    sky_tbl = Table()
    skyfile =f'/data1/Skycatalog/APASS/{target}.csv'
    sky_table = ascii.read(skyfile)
    if filter_ in ['V','R']:
        try:
            sky_table = SDSS_to_JH(sky_table)
            sky_table = sky_table[sky_table[f'e_{filter_}_mag']<0.1]
            sky_tbl['RAJ2000'] = sky_table['RAJ2000']
            sky_tbl['DEJ2000'] = sky_table['DEJ2000']
            sky_tbl[f'{filter_}_mag'] = sky_table[f'{filter_}_mag']
            sky_tbl[f'e_{filter_}_mag'] = sky_table[f'e_{filter_}_mag']
            print(f'{target} mathing with APASS')
            refcat = 'APASS'
        except:
            print(f'{target} not found in APASS')
            refcat = None
    else:
        try:
            sky_table = sky_table[sky_table[f'e_{filter_}_mag']<0.1]
            sky_tbl['RAJ2000'] = sky_table['RAJ2000']
            sky_tbl['DEJ2000'] = sky_table['DEJ2000']
            sky_tbl[f'{filter_}_mag'] = sky_table[f'{filter_}_mag']
            sky_tbl[f'e_{filter_}_mag'] = sky_table[f'e_{filter_}_mag']
            print(f'{target} mathing with APASS')
            refcat = 'APASS'
        except:
            print(f'{target} not found in APASS')
            refcat = None
    return sky_tbl, refcat

def LOAD_SMSS(target, filter_):
    from conversion import SDSS_to_JH
    from astropy.io import ascii
    from astropy.table import Table
    
    sky_tbl = Table()
    try:
        skyfile =f'/data1/Skycatalog/Skymapper/{target}.csv'
        sky_table = ascii.read(skyfile)
        sky_table =sky_table[
            (sky_table['flags'] == 0)&
            (sky_table['nch_max'] ==1)&
            (sky_table['g_ngood'] >=1)&
            (sky_table['r_ngood'] >=1)&
            (sky_table['i_ngood'] >=1)&
            (sky_table['g_psf'] <18)&
            (sky_table['r_psf'] < 18)&
            (sky_table['i_psf'] < 17)&
            (sky_table['g_psf'] > 13)&
            (sky_table['r_psf'] > 13)&
            (sky_table['i_psf'] > 11)&
            (sky_table['class_star']>0.90)&
            (sky_table['e_g_psf'] <0.1)&
            (sky_table['e_r_psf'] <0.1)&
            (sky_table['e_i_psf'] <0.1)
            ]
        sky_tbl['RAJ2000'] = sky_table['raj2000']
        sky_tbl['DEJ2000'] = sky_table['dej2000']
        sky_tbl[f'{filter_}_mag'] = sky_table[f'{filter_}_psf']
        sky_tbl[f'e_{filter_}_mag'] = sky_table[f'e_{filter_}_psf']
        print(f'{target} mathing with SkymapperDR3')
        refcat = 'SMSS'
    except:
        print(f'{target} not found in SMSS')
        refcat = None
    return sky_tbl, refcat

def LOAD_PS1(target, filter_):
    from conversion import PANSTARRS1_to_JH
    from astropy.io import ascii
    from astropy.table import Table
    
    sky_tbl = Table()
    skyfile =f'/data1/Skycatalog/PanSTARRS1/{target}.csv'
    sky_table = ascii.read(skyfile)
    
    if filter_ in ['V','R']:
        try:
            sky_table = PANSTARRS1_to_JH(sky_table)
            sky_table = sky_table[
                (sky_table[f'e_{filter_}_mag']<0.1)&
                (sky_table[f'{filter_}_mag']<18)&
                (sky_table[f'{filter_}_mag']>14)&
                (sky_table[f'{filter_}_mag']-sky_table[f'{filter_}_Kmag']<0.05)
                ]
            sky_tbl['RAJ2000'] = sky_table['RAJ2000']
            sky_tbl['DEJ2000'] = sky_table['DEJ2000']
            sky_tbl[f'{filter_}_mag'] = sky_table[f'{filter_}_mag']
            sky_tbl[f'e_{filter_}_mag'] = sky_table[f'e_{filter_}_mag']
            print(f'{target} mathing with PS1')
            refcat = 'PS1'
        except:
            print(f'{target} not found in PS1')
            refcat = None
    else:
        try:
            sky_table = sky_table[
                (sky_table[f'e_{filter_}mag']<0.1)&
                (sky_table[f'{filter_}mag']<18)&
                (sky_table[f'{filter_}mag']>14)&
                (sky_table[f'{filter_}mag']-sky_table[f'{filter_}Kmag']<0.05)
                ]
            sky_tbl['RAJ2000'] = sky_table['RAJ2000']
            sky_tbl['DEJ2000'] = sky_table['DEJ2000']
            sky_tbl[f'{filter_}_mag'] = sky_table[f'{filter_}mag']
            sky_tbl[f'e_{filter_}_mag'] = sky_table[f'e_{filter_}mag']
            print(f'{target} mathing with PS1')
            refcat = 'PS1'
        except:
            print(f'{target} not found in APASS')
            refcat = None
    return sky_tbl, refcat

def SEXTRACTOR(image, conf_param):
    from astropy.io import fits
    from astropy.io import ascii
    
    hdu = fits.open(image)[0]
    hdr = hdu.header
    target = hdr['OBJECT']


    config = ''
    for param in conf_param.keys():
        config += f'-{param} {conf_param[param]} '

    os.chdir('/data2/sextractor/')
    os.system(f'source-extractor {image} {config}')
    result = ascii.read('zeropoint.cat')
    
    return result

def CROSS_MATCH(obj_catalog,sky_catalog,max_distance_second=5):
    from astropy.coordinates import match_coordinates_sky
    # input = obj_catalog(SkyCoord), sky_catalog(SkyCoord), max_distance_second(match threshold)
    # output = matched_obj_idx, matched_cat_idx, no_matched_obj_idx
    closest_ids, closest_dists, closest_dist3d = match_coordinates_sky(obj_catalog,sky_catalog)
    max_distance = max_distance_second/3600
    matched_object_idx = []
    matched_catalog_idx =[]
    no_matched_object_idx = []
    for i in range(len(closest_dists)):
        if closest_dists.value[i] < max_distance:
            matched_object_idx.append(i)
            matched_catalog_idx.append(closest_ids[i])
        else:
            no_matched_object_idx.append(i)
    return matched_object_idx, matched_catalog_idx,no_matched_object_idx

#%%


def UPDATE_IMAGE(image, obs_info, refcat = 'PS1', class_star_cut = 0.9, magupper = 12, maglower = 16, check = True):
    
    import numpy as np
    import os
    from astropy.io import fits, ascii
    from astropy.coordinates import SkyCoord
    import matplotlib.pyplot as plt
    from astropy.stats import sigma_clip
    from astropy.table import Table
    import glob
            
    #source extractor configuration parameters
    
    conf_param = dict(
        # Default configuration file for Source Extractor 2.25.0
        # EB 2018-02-08
        #
         
        #-------------------------------- Catalog ------------------------------------
         
        CATALOG_NAME     ='test.cat',       # name of the output catalog
        CATALOG_TYPE     ='ASCII_HEAD',     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
                                        # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC
        PARAMETERS_NAME  ='default.param',  # name of the file containing catalog contents
         
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
         
        PHOT_APERTURES   ='6',              # MAG_APER aperture diameter(s) in pixels
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
         
        SEEING_FWHM      =5,            # stellar FWHM in arcsec
        STARNNW_NAME     ='default.nnw',    # Neural-Network_Weight table filename
         
        #------------------------------ Background -----------------------------------
         
        BACK_TYPE        ='AUTO',           # AUTO or MANUAL
        BACK_VALUE       =0.0,            # Default background value in MANUAL mode
        BACK_SIZE        =64,             # Background mesh: <size> or <width>,<height>
        BACK_FILTERSIZE  =3,              # Background filter: <size> or <width>,<height>
         
        BACKPHOTO_TYPE   ='GLOBAL',         # can be GLOBAL or LOCAL
        BACKPHOTO_THICK  =24,             # thickness of the background LOCAL annulus
        BACK_FILTTHRESH  =0.0,            # Threshold above which the background-
                                        # map filter operates
         
        #------------------------------ Check Image ----------------------------------
         
        CHECKIMAGE_TYPE  ='NONE',           # can be NONE, BACKGROUND, BACKGROUND_RMS,
                                        # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                        # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                        # or APERTURES
        CHECKIMAGE_NAME  ='check.fits',     # Filename for the check-image
         
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
    conf_param['CATALOG_NAME'] = 'zeropoint.cat'
    conf_param['PARAMETERS_NAME'] = 'zeropoint.SEparam'
    conf_param['GAIN'] = float(obs_info['gain'])
    conf_param['PIXEL_SCALE'] = float(obs_info['pixelscale'])
    conf_param['SATUR_LEVEL'] = 65536
            
    print(60*'-')
    print(f'Image processing:{os.path.basename(image)}')
    hdr = fits.getheader(image)
    data = fits.getdata(image)
    observatory = obs_info['obs']
    
    target = hdr['OBJECT']
    filter_ = hdr['FILTER']
    
    pixsize = round(float(obs_info['pixelscale']),3)
    
    bkgsize = TARGETSIZE_PIXEL(target, obs_info, observatory)
    #APASS > SMSS > PS1 Reference catalog 

    if refcat == 'APASS':
        sky_tbl, ref = LOAD_APASS(target, filter_)
    if refcat == 'PS1':
        sky_tbl, ref = LOAD_PS1(target, filter_)
    if refcat == 'SMSS':
        sky_tbl, ref = LOAD_SMSS(target, filter_)
    if not refcat in ['APASS','PS1','SMSS']:
        try:
            sky_tbl, ref = LOAD_PS1(target, filter_)
            if ref == None:
                sky_tbl, ref = LOAD_APASS(target, filter_)
            if ref == None:
                sky_tbl, ref = LOAD_SMSS(target, filter_)
            if ref == None:
                print(f'No reference catalog exists for {target}')
        except:
            pass
    if len(sky_tbl) == 0:
        return print(f'No sky reference found for {refcat}')
    
    sky_tbl = sky_tbl[(sky_tbl[f'{filter_}_mag'] > magupper)&
                      (sky_tbl[f'{filter_}_mag'] < maglower)]

    obj_tbl1 = SEXTRACTOR(image, conf_param)
    sorted_tbl1 = obj_tbl1[(obj_tbl1['FLAGS'] == 0)&
            (obj_tbl1['MAGERR_AUTO']<0.1)&
            (3600*obj_tbl1['FWHM_WORLD']>1)
            ]
    sky_tbl
    obj_tbl1
    sky_coord = SkyCoord(sky_tbl['RAJ2000'],sky_tbl['DEJ2000'],unit = 'deg')
    obj_coord = SkyCoord(sorted_tbl1['ALPHA_J2000'],sorted_tbl1['DELTA_J2000'],unit = 'deg')
    
    matched_obj_idx1, matched_sky_idx1, no_matched_obj_idx1 = CROSS_MATCH(obj_coord, sky_coord, 10)


    if len(matched_obj_idx1) < 5:
        print(f'matching failed for {os.path.basename(image)} with {ref}')

    seeing = round(3600*np.median(sigma_clip(sorted_tbl1[matched_obj_idx1]['FWHM_WORLD'],sigma=3,maxiters=1).data),3)
    
    conf_param['PHOT_APERTURES'] = round(3*seeing/pixsize,3)
    conf_param['SEEING_FWHM'] = seeing
    conf_param['BACK_SIZE'] = bkgsize

    # second SExtractor running
    obj_tbl2 = SEXTRACTOR(image, conf_param)
    
    skysig = round(obj_tbl2['THRESHOLD'][0]/1.5,3)
    skyval  = np.median(data[~np.isnan(data)].flatten())
    sorted_tbl2 = obj_tbl2[(obj_tbl2['FLAGS'] == 0)&
            (obj_tbl2['MAGERR_AUTO']<0.1)&
            (obj_tbl2['CLASS_STAR']>class_star_cut)&
            (3600*obj_tbl2['FWHM_WORLD']>1)
            ]
    sky_coord = SkyCoord(sky_tbl['RAJ2000'],sky_tbl['DEJ2000'],unit = 'deg')
    obj_coord = SkyCoord(sorted_tbl2['ALPHA_J2000'],sorted_tbl2['DELTA_J2000'],unit = 'deg')
    matched_obj_idx2, matched_sky_idx2, no_matched_obj_idx2 = CROSS_MATCH(obj_coord, sky_coord, seeing)
    
    zplist = -sorted_tbl2['MAG_APER'][matched_obj_idx2] + sky_tbl[f'{filter_}_mag'][matched_sky_idx2]
    zpclip = sigma_clip(zplist,sigma=3,maxiters=1)
    zperlist = np.sqrt(sorted_tbl2['MAGERR_APER'][matched_obj_idx2]**2 + sky_tbl[f'e_{filter_}_mag'][matched_sky_idx2]**2)
    seeinglist = sorted_tbl2['FWHM_WORLD'][matched_obj_idx2]

    selected_zplist = zplist[~zpclip.mask].data
    selected_zperlist = zperlist[~zpclip.mask].data
    selected_seeinglist = seeinglist[~zpclip.mask].data
    
    zper = round(np.median(selected_zperlist),4)
    zp = round(np.median(selected_zplist),3)
    depth = round(-2.5*np.log10(5*skysig*np.sqrt(np.pi*((1.5*seeing/pixsize)**2))) + zp,3)
    seeing = round(3600*np.median(selected_seeinglist),3)
    
    if check == True:

        plt.figure(figsize = (10,6))
        plt.subplot(2,1,1)
        plt.title(os.path.basename(image))
        plt.ylabel('Seeing[arcsec]')
        plt.scatter(obj_tbl2['MAG_AUTO']+zp, 3600*obj_tbl2['FWHM_WORLD'], c = 'k', alpha = 0.2, label = 'All')
        plt.scatter(sorted_tbl2['MAG_AUTO'][matched_obj_idx2]+zp, 3600*sorted_tbl2['FWHM_WORLD'][matched_obj_idx2], c = 'r', alpha = 0.3, label = 'Matched')
        plt.scatter(sorted_tbl2['MAG_AUTO'][matched_obj_idx2][~zpclip.mask]+zp, 3600*sorted_tbl2['FWHM_WORLD'][matched_obj_idx2][~zpclip.mask], c = 'r', alpha = 1, label = 'Selected')
        plt.ylim(0,10)
        plt.xlim(8,20)
        plt.legend()
        
        plt.subplot(2,1,2)
        plt.xlabel('Mag[AB]')
        plt.ylabel('ZP[AB]')
        plt.errorbar(sorted_tbl2['MAG_APER'][matched_obj_idx2]+zp, zplist, yerr = zperlist , fmt= 'o',  c = 'r',alpha = 0.3, label = f'ZP=clipped')
        plt.errorbar(sorted_tbl2['MAG_APER'][matched_obj_idx2][~zpclip.mask]+zp, selected_zplist, yerr = selected_zperlist , fmt= 'o',  c = 'r', label = f'ZP={zp}[{zper}]')
        plt.text(9,zp-1,f'Depth : {depth}\n Seeing : {seeing}')
        plt.ylim(zp-2,zp+2)
        plt.grid()
        plt.xlim(8,20)
        plt.legend()
        os.chdir(f'{os.path.dirname(image)}')
        plt.savefig(f'{os.path.basename(image).replace(".fits",".check")}.png')
        plt.show()
    
    hdr['SKYSIG'] = (skysig)
    hdr['UL5_4'] = (depth)
    hdr['SEEING'] = (seeing)
    hdr['ZP_4'] = (zp, 'ZP for Mag Aper(Seeing based)[AB]')
    hdr['ZPER_4'] = (zper, 'ZP Error for Mag Aper(Seeing based)[AB]')
    hdr['REFCAT'] = (ref, 'Reference Catalog')
    hdr['SKYVAL'] = (skyval, 'Background median value')
    ########
    
    os.chdir(f'{os.path.dirname(image)}')
    fits.writeto(os.path.basename(image),data,hdr, overwrite=True)
    
    return depth, seeing, zp, zper
