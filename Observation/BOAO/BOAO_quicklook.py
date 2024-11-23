#%%
import numpy as np
from datetime import datetime
import os
import glob
from matplotlib import pyplot as plt
from astropy.table import Table
from astropy.io import ascii, fits
from scipy.ndimage import median_filter
from skimage.feature import peak_local_max
from astropy.stats import sigma_clip, gaussian_fwhm_to_sigma
from numpy.polynomial.chebyshev import chebfit, chebval
from astropy.modeling.models import Gaussian1D
from astropy.modeling.fitting import LevMarLSQFitter
from matplotlib import gridspec
from ccdproc import ImageFileCollection
from astropy.table import vstack, join
import numpy as np
from scipy.interpolate import interp1d
#%% Define the working directory
WD = '/Users/hhchoi1022/Desktop/0112/'
dpi  = 200
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


def divide_calib_files(filelist: list or np.array,
                       object_key = 'object', # object key name in fits header
                       bias_key = 'BIAS', # bias object key name in fits header
                       dark_key = 'DARK', # dark object key name in fits header
                       flat_key = 'FLAT', # flat object key name in fits header
                       remove = False):
    if len(filelist) == 0:
        raise ValueError('Length of the filelist is Zero')
    info_tbl = get_imginfo(filelist)
    bias_tbl = info_tbl[info_tbl[object_key] == bias_key]
    dark_tbl = info_tbl[info_tbl[object_key] == dark_key]
    flat_tbl = info_tbl[info_tbl[object_key] == flat_key]
    biaspathlist = []
    darkpathlist = []
    flatpathlist = []
    for file_ in bias_tbl['file']:
        biaspath = os.path.dirname(file_)+'/BIAS/'
        biaspathlist.append(biaspath + os.path.basename(file_))
        os.makedirs(biaspath, exist_ok = True)
        if remove:
            os.system(f'mv {file_} {biaspath}')
        else:
            os.system(f'cp {file_} {biaspath}')

    for file_ in dark_tbl['file']:
        darkpath = os.path.dirname(file_)+'/DARK/'
        darkpathlist.append(darkpath + os.path.basename(file_))
        os.makedirs(darkpath, exist_ok = True)
        if remove:
            os.system(f'mv {file_} {darkpath}')
        else:
            os.system(f'cp {file_} {darkpath}')
    for file_ in flat_tbl['file']:
        flatpath = os.path.dirname(file_)+'/FLAT/'
        flatpathlist.append(flatpath + os.path.basename(file_))
        os.makedirs(flatpath, exist_ok = True)
        if remove:
            os.system(f'mv {file_} {flatpath}')
        else:
            os.system(f'cp {file_} {flatpath}')
    return biaspathlist, darkpathlist, flatpathlist

#%%
def open_images_ds9(filelist,
                    shell : str = 'zsh'): # If using bash shell, 'bash'
    import subprocess
    names = ""
    if (type(filelist) == str) | (type(filelist) == np.str_):
        names = filelist
    else:
        for file in filelist:
            names += file+" "
    ds9_command = "ds9 "+names+" &"
    sp = subprocess.Popen(["/bin/%s"%shell, "-i", "-c", ds9_command])
    sp.communicate()
    print('Running "'+ds9_command+'" in the terminal...')

#%%
def cut_and_flip(file_ : str,
                 xmin : int = 640,
                 xmax : int = 950):
    # Cutout and flip the spectrum image
    hdu = fits.open(file_)
    data = hdu[0].data
    header = hdu[0].header
    CUT_data = data[:,xmin:xmax]
    Turn = CUT_data.T
    Flip = np.flip(Turn)
    SAVE_cut = os.path.dirname(file_)+'/cut_'+os.path.basename(file_)
    fits.writeto(SAVE_cut, Flip, header = header, overwrite = True)
    return SAVE_cut

#%%
def make_master_bias(biaslist : list or np.array):
    # Make master bias image
    Master_bias = []
    for i in range(0,len(biaslist)):
        hdul = fits.open(biaslist[i])
        bias_data = hdul[0].data
        bias_data = np.array(bias_data).astype('float64')
        Master_bias.append(bias_data)
    MASTER_Bias = np.median(Master_bias,axis=0)
    # Making header part
    bias_header = hdul[0].header  # fetch header from one of bias images
    bias_header['OBJECT'] = 'Bias'
    # - add information about the combine process to header comment
    bias_header['COMMENT'] = f'{biaslist} bias images are median combined on ' \
                            + datetime.now().strftime('%Y-%m-%d %H:%M:%S (KST)')
    SAVE_bias = WD+'/Master_Bias.fits'  # path to save master bias
    fits.writeto(SAVE_bias, MASTER_Bias, header=bias_header, overwrite=True)
    print('Image combined (PATH=%s)'%SAVE_bias)
    return SAVE_bias
    
# %%
def make_master_dark(darklist : list or np.array,
                     mbiaspath : str):
    # Make master dark image
    Mbias = fits.getdata(mbiaspath)
    def get_exptimeset(filelist):
        exptimeset = []
        for i in range(len(filelist)):
            hdul = fits.open(filelist[i])[0]
            header = hdul.header
            exp = header['EXPTIME']  # get exposure time from its header 
            exptimeset.append(exp)
        exptimeset = set(exptimeset)  # only unique elements will be remained
        exptimeset = sorted(exptimeset)
        return exptimeset
    exptimeset = get_exptimeset(filelist = darklist)

    # Making master dark image for each exposure time
    SAVE_darklist = []
    for exp_i in exptimeset:
        Master_dark = []
        for i in range(len(darklist)):
            hdul = fits.open(darklist[i])[0]
            header = hdul.header
            exp = header['EXPTIME']

            if exp == exp_i:
                data = hdul.data
                data = np.array(data).astype('float64')
                bdata = data - Mbias  # bias subtracting
                Master_dark.append(bdata)
                
        MASTER_dark = np.median(Master_dark, axis=0)  # median combine bias-subtracted dark frames
        header['COMMENT'] = 'Bias_subtraction is done' + datetime.now().strftime('%Y-%m-%d %H:%M:%S (KST)')
        header["COMMENT"] = f'{len(Master_dark)} dark images are median combined on '\
                            + datetime.now().strftime('%Y-%m-%d %H:%M:%S (KST)')
        SAVE_dark = WD+'Master_Dark_'+str(exp_i)+'s.fits'
        fits.writeto(SAVE_dark, MASTER_dark, header=header, overwrite=True)
        print('Image combined (PATH=%s)'%SAVE_dark)
        SAVE_darklist.append(SAVE_dark)
    return SAVE_darklist

# %%
def make_master_flat(flatlist : list or np.array,
                     mbiaspath : str,
                     mdarkpath : str = None):

    Mbias = fits.getdata(mbiaspath)
    if isinstance(mdarkpath, str):
        Mdark = fits.getdata(mdarkpath)
        dark_info = fits.getheader(mdarkpath)
        dark_exptime = float(dark_info['EXPTIME'])
    
    # Make master Flat image
    Master_flat = []
    for i in range(len(flatlist)):
        hdul = fits.open(flatlist[i])[0]
        header = hdul.header
        flat_exptime = float(header['EXPTIME'])
        
        data = hdul.data
        data = np.array(data).astype('float64')
        
        bdata = data - Mbias # bias subtraction
        dbdata = bdata 
        header['COMMENT'] = 'Bias_subtraction is done' + datetime.now().strftime('%Y-%m-%d %H:%M:%S (KST)')
        if isinstance(mdarkpath, str):
            dbdata = bdata - Mdark*flat_exptime/dark_exptime  # dark subtraction
            header["COMMENT"] = 'Dark_subtraction is done'+ datetime.now().strftime("%Y-%m-%d %H:%M:%S (KST)") + 'with'+ mdarkpath
        Master_flat.append(dbdata)
        
    MASTER_flat = np.median(Master_flat,axis=0)
    flat = np.mean(MASTER_flat,axis=0)
    sflat = median_filter(flat, 30) # median filtering scale should be defined interactively
    nor_flat2d = []
    for i in range(MASTER_flat.shape[0]):
        flat = MASTER_flat[i,:]
        sflat = median_filter(flat, 30) # this might work better for cutout frames
        nor_flat2d.append(flat / sflat)
    nor_flat2d = np.array(nor_flat2d)    
    header['COMMENT'] = 'Normalized' + datetime.now().strftime('%Y-%m-%d %H:%M:%S (KST)')
    SAVE_flat = WD+'Master_Flat.fits'
    fits.writeto(SAVE_flat, nor_flat2d, header=header, overwrite=True)
    print('Image combined (PATH=%s)'%SAVE_flat)
    return SAVE_flat
    
# %%
def preprocess(file_ : str,
               mbiaspath : str,
               mflatpath : str,
               mdarkpath : str = None):
    
    Mbias = fits.getdata(mbiaspath)
    Mflat = fits.getdata(mflatpath)
    if isinstance(mdarkpath, str):
        Mdark = fits.getdata(mdarkpath)
        dark_info = fits.getheader(mdarkpath)
        dark_exptime = dark_info['EXPTIME']
    
    hdul = fits.open(file_)[0]
    data = hdul.data
    header = hdul.header
    object_exptime = header['EXPTIME']
    bdata = data - Mbias
    header['COMMENT'] = 'Bias_subtraction is done' + datetime.now().strftime('%Y-%m-%d %H:%M:%S (KST)')
    if isinstance(mdarkpath, str):
        dbdata = bdata - Mdark*float(object_exptime)/float(dark_exptime)
        header["COMMENT"] = 'Dark_subtraction is done'+ datetime.now().strftime("%Y-%m-%d %H:%M:%S (KST)") + 'with'+ mdarkpath
    else:
        dbdata = bdata
    fdbjdata = dbdata/Mflat
    header['COMMENT'] = 'Flat correction is done'+ datetime.now().strftime("%Y-%m-%d %H:%M:%S (KST)") + 'with'+ mflatpath
    SAVE_object = WD + "Calib_" + os.path.basename(file_)
    fits.writeto(SAVE_object, fdbjdata, header = header, overwrite = True)
    print('Done! : %s'%SAVE_object)
    return SAVE_object

#%%
def get_spec_source(spec : np.array,
                    dispersion_stepsize : int = 10,
                
                    do_guess_peak : bool = True,
                    sky_lower_min_from_peak = 20,
                    sky_lower_max_from_peak = 40,
                    sky_upper_min_from_peak = 20,
                    sky_upper_max_from_peak = 40,
                    guess_peak : int = 147,
                    guess_fwhm : float = 5,
                    peak_search_width : int = 15,
                    spectrometry_aperture_size : int = 3,
                    visualize = True,
                    save_fig = False,
                    save_name : str = '',
                    do_spectrometry : bool = True,
                    
                    fit_aperture_order : int = 9,
                    fit_aperture_sigma : int = 3,
                    fit_aperture_Niter : int = 10,
                    fit_sky_order : int = 2,
                    fit_sky_sigma : int = 3,
                    fit_sky_Niter : int = 5
                    ) -> np.array:
    '''
    -----
    Usage : Detect the array from a 2D spectrum image(X axis must be dispersion axis)
            & Get the spectrum of the source by aperture tracing with Chebyshev polynomial function
    
    When there are "multiple sources" in the spectrum image:
    set 
        1. do_guess_peak = False
        2. guess_peak = manual pixel coordinate of the source(spatial axis) 
        3. peak_search_width = manual value sufficient to exclude uninterested sources
    otherwise(when only one source exists):
    you can manually set the source center with the above process 
    or 
    just set the do_guess_peak = True
    -----
    Parameters
    1. spec : np.array
        2D array of the spectrum image
    2. dispersion_stepsize : int = 10
        Stepsize along the dispersion axis for aperture tracing
    3. do_guess_peak : bool = True
        If True, guess the source peak coordinate automatically. Do not use when multiple sources exist
    4. sky_lower_min_from_peak : int = 50
        Determine the lower start coordinate of the sky from the peak
    5. sky_lower_max_from_peak : int = 80
        Determine the lower end coordinate of the sky from the peak
    6. sky_upper_min_from_peak : int = 50
        Determine the upper start coordinate of the sky from the peak
    7. sky_upper_max_from_peak : int = 80
        Determine the upper end coordinate of the sky from the peak
    8. guess_peak : int = 150
        (Only when do_guess_peak == False) Manually set the guess of the peak coordinate of the source
    9. guess_fwhm : float = 5
        Set the guess of the fwhm of the source in pixel
    10. peak_search_width : int = 15
        Set the searching area of the source. With guess_fwhm, searching area will be [guess_peak-peak_search_width, guess_peak+peak_search_width]
    11. spectrometry_aperture_size : int = 3
        Set the aperture size of the source for spectrometry in sigma. For example, when setting as 3, 3*seeing aperture is set
    12. visualize : bool = True
        If True, visualize the process for spectrometry(aperture tracing and result aperture)
    13. do_spectrometry : bool = True
        If True, get spectrum
    14. (Advanced)fit_aperture_order : int = 9
        Set the order of Chebyshev polynomial when aperture tracing
    15. (Advanced)fit_aperture_sigma : int = 3
        Set the sigma value for sigma clipping when aperture tracing 
    16. (Advanced)fit_aperture_Niter : int = 5
        Set the number of maximum iteraction for Chebyshev fitting when aperture tracing  
    17. (Advanced)fit_sky_order : int = 9
        Set the order of Chebyshev polynomial when sky fitting
    18. (Advanced)fit_sky_sigma : int = 3
        Set the sigma value for sigma clipping when sky fitting
    19. (Advanced)fit_sky_Niter : int = 5
        Set the number of maximum iteraction for Chebyshev fitting when sky fitting
    -----
    Return
    1. (Only when do_spectrometry == True) 1D spectrum of the source 
    2. SNR of the source : np.array 
        SNR is calculated with the sigma value of the sky region and the intensity of the source
        (WARNING) This is not SNR of the specific Line!
    '''
    
    N_AP = len(spec[0])//dispersion_stepsize
    fitter = LevMarLSQFitter()
    if do_guess_peak:
        apall_1 = np.sum(spec[:,len(spec[0])//2-dispersion_stepsize:len(spec[0])//2+2*dispersion_stepsize],axis=1)
        peak_pix = peak_local_max(apall_1, num_peaks=10,
                                min_distance = 30,
                                threshold_abs=np.median(apall_1))
        peak = []
        for i in peak_pix:
            peak.append(apall_1[i])
        order = peak.index(max(peak))
        pix_peak = peak_pix[order] 
        guess_peak = pix_peak[0] # center
    
    # Define the sky
    ap_sky = np.array([guess_peak-sky_lower_max_from_peak,guess_peak-sky_lower_min_from_peak,
                   guess_peak+sky_upper_min_from_peak,guess_peak+sky_upper_max_from_peak]) # setting the sky area
    x_sky = np.hstack((np.arange(ap_sky[0], ap_sky[1]), 
                    np.arange(ap_sky[2], ap_sky[3])))  # bring the 
    aptrace = []
    aptrace_fwhm = []
    for i in range(N_AP - 1):
        lower_cut, upper_cut = i*dispersion_stepsize, (i+1)*dispersion_stepsize
        apall_i = np.sum(spec[:, lower_cut:upper_cut], axis=1)
        x = np.arange(0, len(apall_i))
        sky_val = np.hstack( (apall_i[ap_sky[0]:ap_sky[1]], 
                        apall_i[ap_sky[2]:ap_sky[3]]) )
        clip_mask = sigma_clip(sky_val,
                            sigma=3,
                            maxiters=5).mask
        coeff, _ = chebfit(x_sky[~clip_mask], 
                           sky_val[~clip_mask],
                           deg=fit_sky_order,
                           full=True)
        apall_i -= chebval(x,coeff)  # Object profile - the fitted sky

    
        search_min = int(guess_peak - peak_search_width)
        search_max = int(guess_peak + peak_search_width)
        cropped = apall_i[search_min:search_max]
        x_cropped = np.arange(len(cropped))
        
        peak_pix = peak_local_max(cropped,
                                min_distance=peak_search_width - 3,
                                num_peaks=1)

        if len(peak_pix)==0: # return NaN (Not a Number) if there is no peak found. 
            aptrace.append(np.nan)
            aptrace_fwhm.append(0)
            continue
            
        else:
            peak_pix = peak_pix[0][0]   
            g_init = Gaussian1D(amplitude=cropped[peak_pix], # Gaussian fitting to find centers
                                mean=peak_pix,
                                stddev=guess_fwhm * gaussian_fwhm_to_sigma,
                                bounds={'amplitude':(0, 2*cropped[peak_pix]) ,
                                        'mean':(peak_pix-peak_search_width, peak_pix+peak_search_width),
                                        'stddev':(0.00001, 4*guess_fwhm*gaussian_fwhm_to_sigma)})
            fitted = fitter(g_init, x_cropped, cropped)
            center_pix_new = fitted.mean.value + search_min
            aptrace_fwhm.append(fitted.fwhm)
            aptrace.append(center_pix_new)  
    aptrace = np.array(aptrace)
    aptrace_fwhm = np.array(aptrace_fwhm)   
     
    if visualize:
        # Plot the center of profile peak
        plt.figure(figsize=(10,5), dpi = dpi)
        plt.title('Source detection')
        plt.imshow(spec,vmin=np.mean(spec)-2*np.std(spec),vmax=np.mean(spec)+2*np.std(spec))
        plt.plot(np.arange(len(aptrace))*dispersion_stepsize, aptrace,ls='',linewidth = 0.1, marker='x', ms=1,color='r', alpha = 0.5, label = 'Detected')
        plt.xlabel('Dispersion axis',fontsize=15)
        plt.ylabel('Spatial axis',fontsize=15)
        plt.legend()
        if save_fig:
            plt.savefig( WD + f'{save_name}_detection.png')
        #plt.show()

        
    if do_spectrometry:
        nan_idx = np.isnan(aptrace)
        x_aptrace = np.arange(N_AP-1) * dispersion_stepsize
        coeff_aptrace = chebfit(x_aptrace[~nan_idx], aptrace[~nan_idx], deg=fit_aperture_order)
        resid_mask = sigma_clip(aptrace - chebval(x_aptrace, coeff_aptrace), 
                        sigma=fit_aperture_sigma, maxiters=fit_aperture_Niter).mask
        x_aptrace_fin = x_aptrace[~resid_mask]
        aptrace_fin = aptrace[~resid_mask]
        coeff_aptrace_fin = chebfit(x_aptrace_fin, aptrace_fin, deg=fit_aperture_order)   

        fit_aptrace_fin   = chebval(x_aptrace_fin, coeff_aptrace_fin)
        resid_aptrace_fin = aptrace_fin - fit_aptrace_fin
        del_aptrace = ~np.in1d(x_aptrace, x_aptrace_fin) # deleted points #x_aptrace에서 x_aptrace_fin이 없으면 True
        if visualize:
            # Plot the Fitted line & residual
            plt.figure(figsize=(10,8))
            plt.title('Aperture trace')
            gs = gridspec.GridSpec(3, 1)
            ax1 = plt.subplot(gs[0:2])
            ax2 = plt.subplot(gs[2])

            ax1.plot(x_aptrace, aptrace, ls='', marker='+', ms=10,color='lightskyblue', label = 'Detected')
            ax1.plot(x_aptrace_fin, fit_aptrace_fin, ls='--',color='crimson',zorder=10,lw=2,
                    label="Aperture Trace function ({:d}/{:d} used)".format(len(aptrace_fin), N_AP-1))
            ax1.plot(x_aptrace[del_aptrace], aptrace[del_aptrace], ls='', marker='x',color='salmon', ms=10)
            ax1.set_ylabel('Found object position')
            ax1.grid(ls=':')
            ax1.legend()

            ax2.plot(x_aptrace_fin, resid_aptrace_fin, ls='', marker='+')
            ax2.axhline(+np.std(resid_aptrace_fin, ddof=1), ls=':', color='k')
            ax2.axhline(-np.std(resid_aptrace_fin, ddof=1), ls=':', color='k', 
                        label='residual std')
            ax2.set_ylabel('Residual (pixel)')
            ax2.set_xlabel('Dispersion axis (pixel)')
            ax2.grid(ls=':')
            ax2.set_ylim(-3, 3)
            ax2.legend()
            if save_fig:
                plt.savefig( WD + f'{save_name}_aptrace.png')
            #plt.show()
            
            
        ap_fwhm = np.median(aptrace_fwhm[~resid_mask]) # [pix]
        ap_sigma = ap_fwhm *gaussian_fwhm_to_sigma # [pixel/sigma]
        x_ap = np.arange(len(spec[0])) # Pixel along the dispersion axis
        y_ap = chebval(x_ap, coeff_aptrace_fin) # Center of peak for each line
        ap_summed  = []
        ap_err = []
        ap_snr = []
        for i in range(len(spec[0])):
            cut_i = spec[:,i] # Cut spatial direction
            peak_i = y_ap[i]
            ap_sky_i = np.array([int(peak_i)-sky_lower_max_from_peak,int(peak_i)-sky_lower_min_from_peak,
                   int(peak_i)+sky_upper_min_from_peak,int(peak_i)+sky_upper_max_from_peak]) # setting the sky area
            
            # aperture size = apsum_sigma_lower * ap_sigma
            x_obj_lower = int(np.around(peak_i - (2 * ap_sigma)//2)) 
            x_obj_upper = int(np.around(peak_i + (2 * ap_sigma)//2))         
            x_obj = np.arange(x_obj_lower, x_obj_upper)
            obj_i = cut_i[x_obj_lower:x_obj_upper]
            
            # Fitting Sky value
            x_sky = np.hstack( (np.arange(ap_sky_i[0], ap_sky_i[1]),
                                np.arange(ap_sky_i[2], ap_sky_i[3])) )
            sky_val = np.hstack( (cut_i[ap_sky_i[0]:ap_sky_i[1]],
                                cut_i[ap_sky_i[2]:ap_sky_i[3]]) )
            clip_mask = sigma_clip(sky_val, sigma=fit_sky_sigma,
                                maxiters=fit_sky_Niter).mask
            coeff = chebfit(x_sky[~clip_mask],
                            sky_val[~clip_mask],
                            deg=fit_sky_order)

            # Subtract the sky
            sky_i = sky_val - chebval(x_sky, coeff)
            sky_sigma = np.std(sky_i)
            sub_obj_i = obj_i - chebval(x_obj, coeff) # obj - lsky  subtraction
            obj_sigma = np.std(np.mean(sub_obj_i))
            aper_sigma = np.sqrt(sky_sigma**2+obj_sigma**2)
            summed = np.sum(sub_obj_i)
            SNR = summed/aper_sigma
            # Calculate error
            #sig_i = RN **2 + sub_obj_i + chebval(x_obj,coeff)
            # RN**2 + flux_i + sky value 
            
            ap_summed.append(np.sum(sub_obj_i))
            ap_err.append(aper_sigma)
            ap_snr.append(SNR)
            #ap_sig.append(np.sqrt(np.sum(sig_i)))
            
        ap_summed = np.array(ap_summed)
        if visualize:
            # Plot the center of profile peak
            plt.figure(figsize=(10,5), dpi = dpi)
            plt.title('Source aperture')
            plt.imshow(spec,vmin=np.mean(spec)-2*np.std(spec),vmax=np.mean(spec)+2*np.std(spec))
            plt.fill_between(x_ap, y_ap+sky_upper_min_from_peak, y2 = y_ap+sky_upper_max_from_peak, alpha = 0.3, color= 'b')
            plt.fill_between(x_ap, y_ap-sky_lower_max_from_peak, y2 = y_ap-sky_lower_min_from_peak, alpha = 0.3, color= 'b', label  = 'Sky')
            plt.plot(x_ap, y_ap,ls='',linewidth = 0.1, marker='x', ms=1,color='r', alpha = 0.1)
            plt.fill_between(x_ap, y_ap+spectrometry_aperture_size//2 * ap_sigma, y2 = y_ap-spectrometry_aperture_size//2 * ap_sigma, alpha = 0.3, color= 'r', label = '[%d*SEEING]Aperture[%.1fpix]'%(spectrometry_aperture_size, ap_sigma * spectrometry_aperture_size))
            plt.xlabel('Dispersion axis',fontsize=15)
            plt.ylabel('Spatial axis',fontsize=15)
            plt.legend()
            if save_fig:
                plt.savefig( WD + f'{save_name}_aperture_sky.png')
            #plt.show()
        return ap_summed, ap_err, ap_snr

#%%
def err_weighted_combine(all_wave, all_flux, all_error, dispersion=None,
                         **kwargs):
    """
    Weigthed Combine
    - from costool procedure
      COADDITION: flux elements are coadded pixel by pixel according to several
      methods set by the keyword method. In our extensive testing, modified
      exposure time weighting with flagging seems to produce the best
      coadditions with the least spurious features. Thus we have made method=1
      the default.
        method = -1 - simple mean of pixel values, gives too much weight to
                      exposures
        method =  1 - modified exposure weighting: exposure time is modified at
                      each pixel location by flanging and wire shadows
                      (if selected).
        method =  2 - err^-2 weighting, allows error array tuning, but tends to
                      weight toward lower-flux pixels.
        method =  3 - (S/N)^2 weighting, allows for error array tuning, but
                      tends to weight toward higher-flux pixels.

    This code corresponds to method=2, the error weighted combine.
    
    This code is based on the repo below.
    https://github.com/hamogu/spectrum/blob/master/spectrum/coadd.py

    Parameters
    ----------
    all_wave : numpy array
        stacked wavelength array.
    all_flux : numpy array
        stacked flux array.
    all_error : numpy array
        stacked error array.
    dispersion : TYPE, optional
        dispersion (wavelength array) for result spectra. The default is None,
        which selects first wavelength array of given spectra.
    
    **kwargs : kwargs for scipy.interpolation.interp1d

    the number of spectra to be combined should be identical for all_wave,
    all_flux, all_error.

    Returns
    -------
    error-weighted spectrum (wavelength, flux, error)

    """
    n_spec = len(all_flux)
    
    if dispersion is None:
        dispersion = all_wave[0]
    
    spec_shape = (n_spec,len(dispersion))
    fluxes, errors = np.ma.zeros(spec_shape), np.ma.zeros(spec_shape)
    
    for i in range(n_spec):
        f_interp = interp1d(all_wave[i], all_flux[i], **kwargs)
        e_interp = interp1d(all_wave[i], all_error[i], **kwargs)
        f_new, e_new = f_interp(dispersion), e_interp(dispersion)
        fluxes[i,:] = f_new
        errors[i,:] = e_new
        
    # First, make sure there is no flux defined if there is no error.
    errors = np.ma.fix_invalid(errors)
    if np.ma.is_masked(errors):
        fluxes[errors.mask] = np.ma.masked
    # This can be simplified considerably as soon as masked quantities exist.
    fluxes = np.ma.fix_invalid(fluxes)
    # There are no masked quantities yet, so make sure they are filled here.
    flux = np.ma.average(fluxes, axis=0, weights = 1./errors**2.).filled(np.nan)
    error = np.sqrt(1. / np.ma.sum(1./errors**2., axis=0).filled(np.nan))
    
    return dispersion, flux, error
#%% 0. Setting the range of the spectrum & Cutout #######################
cut_xmin = 720    ####### CHANGE HERE
cut_xmax = 1020    ####### CHANGE HERE
raw_filelist = glob.glob(WD+'L*.fits')
cut_filelist = []
for raw_file in raw_filelist:
    cut_filelist.append(cut_and_flip(raw_file, xmin = cut_xmin, xmax = cut_xmax))
    
#%% 0. Make calibration images
Bias_filelist, Dark_filelist, Flat_filelist = divide_calib_files(cut_filelist)
# Master Bias
mbias_path = make_master_bias(Bias_filelist)
# Master Dark
mdark_path = make_master_dark(Dark_filelist, mbiaspath  =mbias_path)
# Master Flat 
if isinstance(mdark_path, str):
    mflat_path = make_master_flat(Flat_filelist, mbiaspath = mbias_path)
else:
    mflat_path = make_master_flat(Flat_filelist, mbiaspath = mbias_path)

#%% 1. Single image spectrum extraction #######################
obj_file = '/Users/hhchoi1022/Desktop/0112/LS20230112_000486.fits'
processed_filelist = []
cut_file = cut_and_flip(obj_file, xmin = cut_xmin, xmax = cut_xmax)
processed_file = preprocess(cut_file, mbiaspath = mbias_path, mflatpath  =mflat_path, mdarkpath = mdark_path[0])
open_images_ds9([mbias_path, mflat_path, mdark_path[0], cut_file, processed_file])

obj = fits.getdata(processed_file)
spec, specerr, snr = get_spec_source(spec = obj, 
                                     sky_lower_min_from_peak = 20,
                                     sky_lower_max_from_peak = 150,
                                     sky_upper_min_from_peak = 20,
                                     sky_upper_max_from_peak = 50,
                                     fit_sky_order=2,
                                     visualize = True,  guess_peak = 80, do_guess_peak= True, peak_search_width= 10, spectrometry_aperture_size= 3)
median_filtersize = 1
spec_filtered = median_filter(spec, median_filtersize)
snr_filtered = median_filter(snr, median_filtersize)
plt.figure(figsize = (10,8), dpi  = dpi)
plt.subplots_adjust(left=0.125, bottom=0.1,  right=0.9, top=0.9, wspace=0.2, hspace=0)
ax1 = plt.subplot(211)
ax1.set_ylim(-500, 7000)
ax1.plot(spec_filtered, c= 'k')
ax1.set_ylabel('Total counts',fontsize=15)
plt.xticks(visible=False)
ax2 = plt.subplot(212, sharex = ax1)
ax2.plot(snr_filtered, c= 'r')
ax2.set_ylim(-100, 300)
ax2.set_xlabel('Dispersion axis',fontsize=15)
ax2.set_ylabel('SNR',fontsize=15)


#%% 2. Multiple images (Combine process included)  spectrum extraction #######################

obj_file_1 = '/Users/hhchoi1022/Desktop/0112/LS20230112_000486.fits'
#obj_file_2 = '/Users/hhchoi1022/Desktop/0112/LS20230112_000487.fits'
obj_file_3 = '/Users/hhchoi1022/Desktop/0112/LS20230112_000488.fits'
obj_file_4 = '/Users/hhchoi1022/Desktop/0112/LS20230112_000489.fits'
obj_file_5 = '/Users/hhchoi1022/Desktop/0112/LS20230112_000490.fits'
obj_file_6 = '/Users/hhchoi1022/Desktop/0112/LS20230112_000491.fits'
obj_filelist = [obj_file_1, obj_file_3, obj_file_4, obj_file_5, obj_file_6]
cut_filelist = [cut_and_flip(obj_file, xmin =cut_xmin, xmax = cut_xmax) for obj_file in obj_filelist]

open_images_ds9([mbias_path, mflat_path, mdark_path[0]])
processed_filelist = []
for obj_file in obj_filelist:
    cut_file = cut_and_flip(obj_file, xmin = cut_xmin, xmax = cut_xmax)
    processed_file = preprocess(cut_file, mbiaspath = mbias_path, mflatpath  =mflat_path, mdarkpath= mdark_path[0])
    processed_filelist.append(processed_file)
show_filelist = processed_filelist 

speclist = []
specerrlist = []
snrlist = []
wavelist = []
for processed_file in processed_filelist:
    obj = fits.getdata(processed_file)
    #spec, snr = get_spec_source(spec = obj, visualize = True,  do_guess_peak= True, peak_search_width= 10)
    spec, specerr, snr = get_spec_source(spec = obj, sky_lower_min_from_peak = 30, sky_lower_max_from_peak = 80, sky_upper_min_from_peak = 20, sky_upper_max_from_peak = 40)
    wave = np.arange(1, len(spec)+1)
    speclist.append(spec)
    specerrlist.append(specerr)
    wavelist.append(wave)

_, combined_spec, combined_err = err_weighted_combine(np.array(wavelist), np.array(speclist), np.array(specerrlist))
combined_snr = combined_spec/combined_err

median_filtersize = 1
spec_filtered = median_filter(spec, median_filtersize)
snr_filtered = median_filter(snr, median_filtersize)
combined_spec_filtered = median_filter(combined_spec, median_filtersize)
combined_snr_filtered = median_filter(combined_snr, median_filtersize)
plt.figure(figsize = (10,8), dpi  = dpi)
plt.subplots_adjust(left=0.125, bottom=0.1,  right=0.9, top=0.9, wspace=0.2, hspace=0)
ax1 = plt.subplot(211)
ax1.set_ylim(-200, 2000)
#ax1.plot(spec_filtered, c= 'k', label = 'Single')
ax1.plot(combined_spec_filtered, c= 'k', label = 'Combined')
ax1.legend()
ax1.set_ylabel('Total counts',fontsize=15)
plt.xticks(visible=False)
ax2 = plt.subplot(212, sharex = ax1)
ax2.set_ylim(0, 150)
#ax2.plot(snr_filtered, c= 'r', label = 'Single')
ax2.plot(combined_snr_filtered, c= 'r', label = 'Combined')
ax2.set_xlabel('Dispersion axis',fontsize=15)
ax2.set_ylabel('SNR',fontsize=15)
ax2.legend()
# %%

#%% 3. Nodding spectrum extraction ##############################################
image_A_set = ['/Users/hhchoi1022/Desktop/0112/LS20230112_000486.fits', '/Users/hhchoi1022/Desktop/0112/LS20230112_000489.fits', '/Users/hhchoi1022/Desktop/0112/LS20230112_000490.fits'] 
image_B_set = ['/Users/hhchoi1022/Desktop/0112/LS20230112_000487.fits', '/Users/hhchoi1022/Desktop/0112/LS20230112_000488.fits', '/Users/hhchoi1022/Desktop/0112/LS20230112_000491.fits'] 
def nodding(image_A, 
            image_B,
            **kwargs
            ) -> np.array:
    cut_image_A = cut_and_flip(image_A, xmin = cut_xmin, xmax = cut_xmax)
    cut_image_B = cut_and_flip(image_B, xmin = cut_xmin, xmax = cut_xmax)
    pre_image_A = preprocess(cut_image_A, mbiaspath = mbias_path, mflatpath= mflat_path)
    pre_image_B = preprocess(cut_image_B, mbiaspath = mbias_path, mflatpath= mflat_path)
    data_A = fits.getdata(pre_image_A)
    data_B = fits.getdata(pre_image_B)
    data_A = data_A.astype('float32')
    data_B = data_B.astype('float32')
    sub_data = data_A - data_B
    spec, specerr, snr = get_spec_source(sub_data, **kwargs)
    return spec, specerr, snr

speclist = []
specerrlist = []
snrlist = []
wavelist = [] 

for i, images in enumerate(zip(image_A_set, image_B_set)):
    image_A, image_B = images
    spec, specerr, snr = nodding(image_A, image_B, sky_upper_min_from_peak = 20, sky_upper_max_from_peak = 40, sky_lower_min_from_peak = 20, sky_lower_max_from_peak = 100)
    wave = np.arange(1, len(spec)+1)
    speclist.append(spec)
    specerrlist.append(specerr)
    snrlist.append(snr)
    wavelist.append(wave)
    
_, combined_spec, combined_err = err_weighted_combine(np.array(wavelist), np.array(speclist), np.array(specerrlist))
combined_snr = combined_spec/combined_err

median_filtersize = 1
spec_filtered = median_filter(combined_spec, median_filtersize)
snr_filtered = median_filter(combined_snr, median_filtersize)
plt.figure(figsize = (10,8), dpi  = dpi)
plt.subplots_adjust(left=0.125, bottom=0.1,  right=0.9, top=0.9, wspace=0.2, hspace=0)
ax1 = plt.subplot(211)
ax1.set_ylim(-500, 7000)
ax1.plot(combined_spec, c= 'k')
ax1.set_ylabel('Total counts',fontsize=15)
#ax1.set_xlim(1800,1850)
plt.xticks(visible=False)
ax2 = plt.subplot(212, sharex = ax1)
ax2.set_ylim(-100, 300)
ax2.plot(combined_spec/combined_err, c= 'r')
ax2.set_xlabel('Dispersion axis',fontsize=15)
ax2.set_ylabel('SNR',fontsize=15)
