

#%%
from astropy.io import fits
from astropy.wcs import WCS
# %%
import glob
imlist = glob.glob('/home/hhchoi1022/Desktop/GRB221009A/KCT_STX16803/Calibration_img/*.fits')
# %%
import os
import numpy as np
def bin_img(image,
            x_bin_factor = 2,
            y_bin_factor = 2,
            prefix = 'bin_',
            method = 'max'):
    hdu = fits.open(image)[0]
    origin_shape = hdu.data.shape
    cut_shape_x = (origin_shape[0] // x_bin_factor) * x_bin_factor
    cut_shape_y = (origin_shape[1] // y_bin_factor) * y_bin_factor
    image_cut = hdu.data[origin_shape[0]-cut_shape_x:, origin_shape[1]-cut_shape_y:]
    image_reshaped = image_cut.reshape(origin_shape[0] // x_bin_factor, x_bin_factor, origin_shape[1] // y_bin_factor ,y_bin_factor)
    if method =='max':
        binned_data = np.max(np.max(image_reshaped, axis= 1), axis = 2)
    if method =='mean':
        binned_data = np.mean(np.mean(image_reshaped, axis= 1), axis = 2)
    if method =='median':
        binned_data = np.median(np.median(image_reshaped, axis= 1), axis = 2)
    binned_data
    fits.writeto((f'{os.path.dirname(image)}/{prefix}{os.path.basename(image)}'), binned_data, hdu.header, overwrite = True)
# %%
imlist
# %%
for image in imlist:
    bin_img(image, x_bin_factor= 3, y_bin_factor= 3, method = 'max')
# %%
