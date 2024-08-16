#%%
from astropy import wcs
from ccdproc import CCDData
from astropy.io import fits
import os
from shutil import ExecError
import astroalign as aa


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


#%%
import glob
imlist = glob.glob('/home/hhchoi1022/Desktop/GRB221009A/KCT_STX16803/i/221013/Calib*.fit')
# %%
refim = imlist[15]
# %%
for image in imlist:
    align_img(image, refim)
# %%
