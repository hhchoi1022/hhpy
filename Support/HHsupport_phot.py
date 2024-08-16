#%%
from shutil import ExecError
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u
import sys
from ccdproc import CCDData
from astropy.stats import sigma_clipped_stats
import os
from ccdproc import Combiner
import numpy as np
from astropy.io import ascii
from ccdproc import ImageFileCollection
from astropy.table import vstack
from astropy.table import Table
from astropy.io import fits
import astroalign as aa
from astropy.time import Time
from astropy import wcs
from astropy.nddata import Cutout2D
import os
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import join 
from time import sleep
from astropy.table import unique
from astropy.visualization import astropy_mpl_style
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import subprocess

#%%



get_obsinfo_file = '/home/hhchoi1022/Desktop/Gitrepo/Config/CCD.dat'
get_tarinfo_file = '/home/hhchoi1022/Desktop/Gitrepo/Config/alltarget.dat'
#%%

def get_imginfo(filelist,
                keywords = 'default'
                ):
    '''
    parameters
    ----------
    1. filelist : str or list or np.array
            (str) filename of the fits image
            (list or np.array) filelist of the fits images
    2. keywords : 
    
    returns 
    -------
    1. result_tbl : astropy.table
            all information in the header of the filelist
    
    notes 
    -----
    "comment" and "history" are removed in the result_tbl due to too long length of the information 
    -----
    '''
    default_keywords = ['jd','group','filter','exptime','object','telescop','instrume','ra','dec','xbinning','ybinning','imagetyp']
    if keywords == 'default':
        keywords = default_keywords
    else:
        keywords = default_keywords + keywords
    if (type(filelist) == str) | (type(filelist) == np.str_):
        filelist = np.array(filelist, ndmin = 1)
    filedirs = set([os.path.dirname(filename) for filename in filelist])
    match_tbl = Table([filelist],names = ['file'])
    result_tbl = Table()
    for filedir in filedirs:
        files = [os.path.basename(file) for file in filelist if os.path.dirname(file)==filedir]
        iminfo = ImageFileCollection(location = filedir, keywords=keywords, filenames = files).summary
        absfiles = []
        for file in iminfo['file']:
            absfiles.append(filedir + '/'+ file)
        iminfo['file'] = absfiles
        result_tbl = vstack([iminfo, result_tbl])
    result_tbl = join(match_tbl, result_tbl, keys = 'file', join_type = 'inner')
    return result_tbl


def get_tarinfo(target,
                key_target = 'obj'):
    
    '''
    parameters
    ----------
    1. target : str
            target name to search in target information file
    2. key_target : str
            the keyword for searching target in target information file
            
    returns 
    -------
    1. target_info : astropy.table
            target information 
    
    notes 
    -----
    This code finds the target in IMSNG targetlist. 
    Output includes basic information of the target
    -----
    '''
    
    all_tarinfo = ascii.read(get_tarinfo_file, format = 'fixed_width')
    if not target in all_tarinfo[key_target]:
        raise AttributeError(f'{target} is not found in target information file')
    else:
        targetinfo = all_tarinfo[all_tarinfo[key_target] == target]
        return targetinfo

def get_obsinfo(observatory, ccd = None, rasamode = 'High',
                key_observatory = 'obs',
                key_ccd = 'ccd'):
    
    '''
    parameters
    ----------
    1. observatory : str
            observatory name to search in observaroy information file
    2. ccd : optional, str
            ccd name of the observatroy (None)
    3. rasamode : optional, str
            rasamode for RASA36 ccd [High, Merge] (High) 
    4. key_observatory : str
            the key for searching observatory in observaroy information file (obs)
    6. key_ccd : str
            the key for seraching ccd in observaroy information file (ccd)
    returns 
    -------
    1. obsinfo : astropy.table
            The information of the observatory/ccd 
    
    notes 
    -----
    -----
    '''
    
    output_valid_func = lambda output: len(output) == 1
    all_obsinfo = ascii.read(get_obsinfo_file, format = 'fixed_width')
    obs_info = all_obsinfo[all_obsinfo[key_observatory] == observatory]
    if not observatory in all_obsinfo[key_observatory]:
        raise AttributeError(f'{observatory} information not exist.\n available :{set(all_obsinfo[key_observatory])}')

    #RASA36
    if observatory == 'RASA36':
        if rasamode == None:
            rasamode = input('RASA36 has multiple modes for observation. Please select one. (High/Merge)')
        if rasamode in ['merge','Merge','MERGE']:
            obs_info = obs_info[obs_info['gain'] > 10]
        if rasamode in ['high', 'High','HIGH']:
            obs_info = obs_info[obs_info['gain'] < 10]
        if output_valid_func(obs_info):
            return obs_info
    # Other observatories (Single CCD)
    if output_valid_func(obs_info):
        return obs_info
    # Other observatories (Multiple CCDs)
    else:
        if ccd == None:
            print('Multiple information of the observatory is found.')
            ccd = input(f'Choose the CCD : {set(obs_info[key_ccd])}')
        obs_info = all_obsinfo[(all_obsinfo[key_observatory] == observatory) & (all_obsinfo[key_ccd] == ccd)]
        if output_valid_func(obs_info):
            return obs_info
    
    if not output_valid_func(obs_info):
        raise AttributeError(f'{ccd} information not exist in {observatory}.\n available CCD name:{list(set(all_obsinfo[all_obsinfo[key_observatory] == observatory][key_ccd]))}')

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

def align_img(target_img,
              reference_img,
              prefix = 'align_'
              ):
    """

    parameters
    ----------
    1. target_img : str
            absolute path of a target image 
    2. reference_img : str
            absolute path of a reference image
            The image as a reference frame. All images will be aligned on the basis of this image
    3. prefix = str
            prefix of the aligned iamge
      
    returns
    ----------
    1. outputname : str
            absolute path of the aligned image
    notes
    ----------
    ----------
    """
    tgt_data, tgt_hdr = fits.getdata(target_img, header = True)
    ref_data, ref_hdr = fits.getdata(reference_img, header = True)
    ref_wcs = wcs.WCS(ref_hdr)
    wcs_hdr = ref_wcs.to_header()
    tgt_hdr.update(wcs_hdr)
    tgt_data = tgt_data.byteswap().newbyteorder()
    ref_data = ref_data.byteswap().newbyteorder()
    
    try:
        aligned_data, footprint = aa.register(tgt_data, ref_data, fill_value = 0, detection_sigma = 3, max_control_points = 30)
        aligned_tgt = CCDData(aligned_data, unit = 'adu', header = tgt_hdr)
        outputname = f'{os.path.dirname(target_img)}/{prefix}{os.path.basename(target_img)}'
        fits.writeto(outputname, aligned_tgt.data, aligned_tgt.header, overwrite = True)
        return outputname
    except:
        raise ExecError('List of matching triangles exhausted before an acceptable transformation was found')

def combine_img(filelist,
                clip = 'sigma',
                combine = 'median',
                scale = 'zero',
                prefix = 'com_',
                refim = None,
                
                # Clipping
                clip_sigma_low = 2,
                clip_sigma_high = 5,
                clip_minmax_min = 3,
                clip_minmax_max = 3,
                clip_extrema_nlow = 1,
                clip_extrema_nhigh = 1,
                
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
    for i, file in enumerate(filelist):
        hdr[f'COMBINE{i+1}'] = os.path.basename(file)
        
    hdr['NCOMBINE'] = int(len(filelist))
    if 'JD' in hdr.keys():
        hdr['JD'] = Time(np.mean([inim.header['JD'] for inim in ccdlist]), format='jd').value
    if 'DATE-OBS' in hdr.keys():
        hdr['DATE-OBS'] = Time(np.mean([Time(inim.header['DATE-OBS']).jd for inim in ccdlist]),format= 'jd').isot

    hdr['EXPTIME'] = int(np.sum([inim.header['EXPTIME'] for inim in ccdlist]))
    
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

    outputname = f'{os.path.dirname(filelist[0])}/{prefix}{os.path.basename(filelist[0])}'
    if (len(filelist) == 1):
            
        ccd.header = hdr
        ccd.write(outputname, overwrite=True, format = 'fits')
    else:
        combined.write(outputname, overwrite=True, format = 'fits')        
    init_mean, init_std = np.mean(ccdlist[0].data), np.std(ccdlist[0].data)
    fin_mean, fin_std = np.mean(combined.data), np.std(combined.data)
    print('Combine complete \n')
    print('Combine information')
    print(60*'=')
    print(f'Ncombine = {len(filelist)}')
    print(f'method   = {clip}(clipping), {combine}(combining)')
    print(f'mean     = {round(init_mean,3)} >>> {round(fin_mean,3)}')
    print(f'std      = {round(init_std,3)} >>> {round(fin_std,3)}')
    print(f'image path = {outputname}')
    
    return outputname

def cutout_img(target_img,
               size = 0.9,
               prefix = 'cut_',
               
               xcenter = None,
               ycenter = None
               ):
    
    '''
    parameters
    ----------
    1. target_img : str
            {abspath} of the target image
    2. size : float or int
            (float) ratio of the cut image (0.9)
            (int) size of the cut image in pixel 
    3. prefix : str
            prefix of the output image
    4. xcenter : optional, int or (float or str)
            (int) pixel coordinate of the center
            (float or str) RA coordinate of the center
    4. ycenter : optional, int or (float or str)
            (int) pixel coordinate of the center
            (float or str) Dec coordinate of the center
    returns 
    -------
    1. outputname : str
            {abspath} of the cutout image
    
    notes 
    -----
    if xcenter, ycenter == None, cutout image will be centered to the center of the original image
    -----
    '''
    
    hdu = fits.open(target_img)[0]
    wcs = WCS(hdu.header)
    if size < 1:
        size = size*int(len(hdu.data))
    if (xcenter == None) & (ycenter == None):
        xcenter, ycenter = len(hdu.data)//2, len(hdu.data)//2
    if (type(xcenter) !=int) & (type(ycenter) !=int):
        center_coords = to_skycoord(xcenter, ycenter)
        cutouted = Cutout2D(data = hdu.data, position = center_coords, size = size, wcs = wcs)
    else:
        cutouted = Cutout2D(data = hdu.data, position = (xcenter, ycenter), size = size, wcs = wcs)
    cutouted_hdu = hdu
    cutouted_hdu.data = cutouted.data
    cutouted_hdu.header.update(cutouted.wcs.to_header())
    cutouted_hdu.header['NAXIS1'] = int(size)
    cutouted_hdu.header['NAXIS2'] = int(size)
    outputname = f'{os.path.dirname(target_img)}/{prefix}{os.path.basename(target_img)}'
    cutouted_hdu.writeto(outputname, overwrite = True)
    return outputname

def subtract_img(target_img,
                 reference_img,
                 prefix = 'sub_',
                 method = 'hotpants',
                 # hotpants config
                 iu = 60000,
                 il = -100000,
                 tu = 600000000,
                 tl = -100000,
                 v = 0,
                 ng = '3 3 1.0 2 0.7 1 0.4'
                ):
    '''
    parameters
    ----------
    1. target_img : str
            {abspath} of the target image to be subtracted
    2. reference_img : str
            {abspath} of the reference image 
    3. prefix : str
            prefix of the output image (sub_)
    4. method : str
            method for subtraction (hotpants)
    5. iu : int
            upper valid data count, image (60000)
    6. tu : int
            upper valid data count, template (600000000)
    7. tl : int
            lower valid data count, template (-100000)
    8. v : int
            level of verbosity, 0-2 (0)
    9. ng : str
            'ngauss degree0 sigma0 .. degreeN sigmaN'
            : ngauss = number of gaussians which compose kernel (3)
            : degree = degree of polynomial associated with gaussian # (3 2 1)
            : sigma  = width of gaussian # (1.0 0.7 0.4)
    returns 
    -------
    
    notes 
    -----
    For more information : https://github.com/acbecker/hotpants
    -----
    '''
    outputname = f'{os.path.dirname(target_img)}/{prefix}{os.path.basename(target_img)}'
    if method == 'hotpants':
        os.system(f'hotpants -c t -n i -inim {target_img} -tmplim {reference_img} -outim {outputname} -iu {iu} -il {il} -tu {tu} -tl {tl} -v {v} -ng {ng} > .out && rm -rf .out')
    return outputname

#%%
def to_skycoord(ra, dec):
    '''
    parameters
    ----------
    1. ra : str or float
            Right ascension in diverse format(see notes)
    2. dec : str or float
            Declination in diverse format(see notes)
    
    returns 
    -------
    1. skycoord : SkyCoord
    
    notes 
    -----
    Current supported formats
        1. 15h32m10s, 50d15m01s
        2. 15 32 10, 50 15 01
        3. 15:32:10, 50:15:01
        4. 230.8875, 50.5369
    -----
    '''
    if (type(ra) == float) | (type(ra) == int) | (type(ra) == str) | (type(ra) == np.str_):
        ra = str(ra) ; dec = str(dec)
        if (':' in ra) & (':' in dec):
            skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg), frame = 'fk5')
        elif ('h' in ra) & ('d' in dec):
            skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg), frame = 'fk5')
        elif (' ' in ra) & (' ' in dec):
            skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg), frame = 'fk5')
        else:
            skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.deg, u.deg), frame = 'fk5')
    else:
        if (type(ra[0]) == str) | (type(ra[0]) == np.str_):
            if (':' in ra[0]) & (':' in dec[0]):
                skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg), frame = 'fk5')
            elif ('h' in ra[0]) & ('d' in dec[0]):
                skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg), frame = 'fk5')
            elif (' ' in ra[0]) & (' ' in dec[0]):
                skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.hourangle, u.deg), frame = 'fk5')
        else:
            skycoord = SkyCoord(ra = ra, dec = dec, unit = (u.deg, u.deg), frame = 'fk5')
    return skycoord

# %%

def run_ds9(filelist):
    '''
    parameters
    ----------
    1. filelist : str or list or np.array
            (str) abspath of the fitsfile for visualization
            (list or np.array) the list of abspath of the fitsfiles for visualization
    
    returns 
    -------
    
    notes 
    -----
    -----
    '''
    ds9_options = "-scalemode zscale -scale lock yes -frame lock image "
    names = ""
    if (type(filelist) == str) | (type(filelist) == np.str_):
        names = filelist
    else:
        for file in filelist:
            names += file+" "
    ds9_command = "ds9 "+ds9_options+names+" &"
    print('Running "'+ds9_command+'" in the terminal...')
    sp = subprocess.Popen(["/bin/zsh", "-i", "-c", ds9_command])
    sp.communicate()
    #os.system(ds9_command)

def select_sources(fitsfile,
                   irafpath = '/data2/iraf'):
    '''
    parameters
    ----------
    1. fitsfile : str
            abspath of the fitsfile for visualization
    2. irafpath : optional, str
            abspath of the iraf folder including iraf configuration file
    
    returns 
    -------
    1. outputpath : str 
            abspath of the output file 
            the output file includes ID, coordinates(ra,dec), note 
    
    notes 
    -----
    This code uses IRAF and ds9 for extracting information and visualization. 
    Default iraf path is set to /data2/iraf
    -----
    '''
    current_dir = os.getcwd()
    try:
        dir_iraf = os.path.abspath(irafpath)    # where 'login.cl' file is locate
        os.chdir(dir_iraf)
        from pyraf import iraf
        os.chdir(current_dir)
    except:
        pass
    
    run_ds9(fitsfile)
    print('============== PRESS "a" to select a source ==============')
    print('============== PRESS "q" to quit ==============')
    sleep(3)
    dir_img = os.path.dirname(fitsfile)
    filename = dir_img+"/tmp.cat"
    if os.path.isfile(filename):
        os.system(f'rm {filename}')
    iraf.imexamine(logfile = filename, keeplog = "yes")
    sources = np.genfromtxt(filename, usecols=(0,1), names=('X','Y'))
    if sources.ndim ==0:
        sources = sources.reshape(1)
    _, hdr = fits.getdata(fitsfile, header = True)
    w = WCS(hdr)
    ralist = []
    declist = []
    result = Table()
    if len(sources) > 0:
        for x, y in sources:
            coords = w.pixel_to_world(x,y)
            ra = coords.ra.value
            dec = coords.dec.value
            ralist.append(ra)
            declist.append(dec)
        IDlist = 1+np.arange(len(ralist))
        result['ID'] = IDlist
        result['ra'] = ralist
        result['dec'] = declist
        result = unique(result, 'ra')
        outputname = dir_img+'/sources.cat'
        result.write(f'{outputname}', format = 'ascii.fixed_width', overwrite = True)
        print(f'============== {len(result)} sources are selected ({len(IDlist)-len(result)} are overlapped))  ==============')
        print(f'============== file path : {outputname}  ==============')
        os.system(f'rm -rf {filename}')
        return outputname
    else:
        return None
    

}

# %%

def printProgress (iteration, total, prefix = '', suffix = '', decimals = 1, barLength = 100): 
    '''
    parameters
    ----------
    1. iteration : int
            iteration counts 
    2. total : int
            total iteration counts
    3. prefix : optional, str 
            prefix for the progress status 
    4. suffix : optional, str 
            suffix for the progress status 
    5. decimals : int
            digit for percentage (1)
    6. barlength : int
            length of the progress bar 
    
    returns 
    -------
    
    notes 
    -------
    This can be used to express current progress for iterative process
    To use this, 
        input {i} as {iteraction} with (i, value in enumerate(value)
        input {len(value)} as {total}
    -----
    '''
    formatStr = "{0:." + str(decimals) + "f}" 
    percent = formatStr.format(100 * (iteration / float(total))) 
    filledLength = int(round(barLength * iteration / float(total))) 
    bar = '#' * filledLength + '-' * (barLength - filledLength) 
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percent, '%', suffix)), 
    if iteration == total: 
        sys.stdout.write('\n') 
    sys.stdout.flush() 
printProgress
# %%

def load_sexconfig():
    '''
    parameters
    ----------

    returns 
    -------
    default_config : dict
            default configuration for Source Extractor
    notes 
    -------
    -----
    '''
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
        DETECT_THRESH    =3,            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
        ANALYSIS_THRESH  =3,            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
         
        FILTER           ='Y',              # apply filter for detection (Y or N)?
        FILTER_NAME      ='default.conv',   # name of the file containing the filter
         
        DEBLEND_NTHRESH  =64,             # Number of deblending sub-thresholds
        DEBLEND_MINCONT  =0.0001,          # Minimum contrast parameter for deblending
         
        CLEAN            ='N',              # Clean spurious detections? (Y or N)?
        CLEAN_PARAM      =1.0,            # Cleaning efficiency
         
        MASK_TYPE        ='CORRECT',        # type of detection MASKing: can be one of
                                        # NONE, BLANK or CORRECT
         
        #-------------------------------- WEIGHTing ----------------------------------
        
        WEIGHT_TYPE      ='NONE',           # type of WEIGHTing: NONE, BACKGROUND,
                                        # MAP_RMS, MAP_VAR or MAP_WEIGHT
        RESCALE_WEIGHTS  ='N',              # Rescale input weights/variances (Y/N)?
        WEIGHT_IMAGE     ='weight.fits',    # weight-map filename
        WEIGHT_GAIN      ='N',              # modulate gain (E/ADU) with weights? (Y/N)
        
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
         
        CHECKIMAGE_TYPE  ='SEGMENTATION',           # can be NONE, BACKGROUND, BACKGROUND_RMS,
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
    return conf_param

def run_sextractor(image, conf_param,
                   sexpath ='/data2/sextractor'):
    '''
    parameters
    ----------
    1. image : str
            abspath of the target image
    2. conf_param : dict
            configuration parameters in dict format. can be load by load_sexconfig()
    3. sexpath : str
            source-extractor directory
    returns 
    -------
    1. result : astropy.table
            source extractor result
    
    notes 
    -------
    -----
    '''
    
    sexconfig = ''
    for param in conf_param.keys():
        sexconfig += f'-{param} {conf_param[param]} '
    currentpath = os.getcwd()
    os.chdir(sexpath)
    os.system(f'source-extractor {image} {sexconfig} >/dev/null 2>&1')
    sexresult = ascii.read(conf_param['CATALOG_NAME'])
    os.chdir(currentpath)
    return sexresult


# %%
def plot_matched(image, 
                obs_ra,
                obs_dec,
                ref_ra,
                ref_dec):
        hdu = fits.open(image)[0]
        wcs = WCS(hdu.header)
        median = np.median(hdu.data[~np.isnan(hdu.data)])
        std = np.std(hdu.data[~np.isnan(hdu.data)])
        vmin = median - 1*std
        vmax = median + 1*std    
        
        
        #plt.style.use(astropy_mpl_style)
        ax = plt.subplot(projection = wcs)
        ax.imshow(hdu.data, origin='lower', vmin = vmin, vmax =vmax, cmap = 'gray')
        obs_ra = obs_ra.degree
        obs_dec = obs_dec.degree
        ax.scatter(obs_ra, obs_dec, transform = ax.get_transform('fk5'), c = 'red', marker = '.',s =10, label = 'Obs')
        ax.scatter(ref_ra, ref_dec, transform = ax.get_transform('fk5'), edgecolor = 'red', facecolor = 'none', s =30, label = 'Ref')
        ax.legend()
# %%
