#%%
from astropy.io import ascii
from astropy.io import fits
import numpy as np
from specutils import Spectrum1D, SpectralRegion
import astropy.units as u
from specutils.manipulation import box_smooth, gaussian_smooth, trapezoid_smooth, extract_region
from specutils.fitting import fit_generic_continuum
from specutils.fitting import estimate_line_parameters
from astropy.modeling import models
from specutils.fitting import fit_lines
from specutils.analysis import equivalent_width, fwzi, gaussian_fwhm, fwhm

spec_b = fits.open('./mbgphH202309180019_uwm.fits')
spec_r = fits.open('./mbgphR202309180019_uwm.fits')
# NIST ASD Output
wl_D1_rest = 5889.95095
wl_D2_rest = 5895.92424
redshift = 0.00424 #NED
# %%
import matplotlib.pyplot as plt
# %%
value_b = spec_b[0].data
header_b = spec_b[0].header

value_r = spec_r[0].data
header_r = spec_r[0].header
value_all = np.append(value_b, value_r[2656:])

# %%
idx_ref_b = int(header_b['CRPIX1'])
wl_ref_b = header_b['CRVAL1']
delwl_ref_b = header_b['CDELT1']

idx_ref_r = int(header_r['CRPIX1'])
wl_ref_r = header_r['CRVAL1']
delwl_ref_r = header_r['CDELT1']

# %%
pix_coord_b = np.arange(idx_ref_b-1, len(value_b))
wl_b = wl_ref_b + pix_coord_b * delwl_ref_b

pix_coord_r = np.arange(0, len(value_r)-idx_ref_r+1)
wl_r = wl_ref_r + pix_coord_r * delwl_ref_r
wl_r = wl_r[-idx_ref_r+1:]

wl_all = np.append(wl_b , wl_r[2656:])

# %%

spec1D_b= Spectrum1D(flux = value_b * u.mJy, spectral_axis = wl_b * u.AA)
spec1D_r= Spectrum1D(flux = value_r * u.mJy, spectral_axis = wl_r * u.AA)
spec1D_all = Spectrum1D(flux = value_all * u.mJy, spectral_axis = wl_all * u.AA)
wl_rest = wl_all / (1+redshift)
spec1D_rest = Spectrum1D(flux = value_all * u.mJy, spectral_axis = wl_rest * u.AA)
# %%
spec1_bsmooth = box_smooth(spec1D_rest, width=10)
spec1_gsmooth = gaussian_smooth(spec1D_rest, stddev=3)
spec1_tsmooth = trapezoid_smooth(spec1D_rest, width=3)
# %% Plot
i = 10000
f, ax = plt.subplots()  
ax.step(spec1_bsmooth.spectral_axis[i:], spec1_bsmooth.flux.value[i:]) 
#ax.set_xlim(5850, 6000)

#%% Slicing
spec_slice = spec1D_rest[5750*u.AA:5960*u.AA]

#%% Smoothing 
spec1_bsmooth = box_smooth(spec_slice, width=30)
f, ax = plt.subplots()  
ax.step(spec1_bsmooth.spectral_axis, spec1_bsmooth.flux.value) 
ax.set_xlim(5750, 5960)
ax.set_ylim(0, 0.003)
#ax.axvline(wl_D1_rest, c ='b', label = 'D1(rest)')
#ax.axvline(wl_D2_rest, c = 'r', label = 'D2(rest)')
plt.legend()
#%% continuum fitting
x = spec1_bsmooth.spectral_axis
g1_fit = fit_generic_continuum(spec1_bsmooth)
y_continuum_fitted = g1_fit(x)
#%% plot
f, ax = plt.subplots()  
ax.step(spec1_bsmooth.spectral_axis, spec1_bsmooth.flux.value) 
ax.step(x, y_continuum_fitted)
# %% Continuum subtraction
spec1D_sub = Spectrum1D(flux = (spec1_bsmooth.flux - y_continuum_fitted), spectral_axis = spec1_bsmooth.spectral_axis)
spec1D_norm = Spectrum1D(flux = (spec1_bsmooth.flux/y_continuum_fitted), spectral_axis = spec1_bsmooth.spectral_axis)
#%% Line finding
from specutils.manipulation import noise_region_uncertainty
noise_region = SpectralRegion(5910*u.AA, 5930*u.AA)
#%%
spectrum = noise_region_uncertainty(spec1D_sub, noise_region)
from specutils.fitting import find_lines_threshold
lines = find_lines_threshold(spectrum, noise_factor=10)  

#%% Extract line regions 
sub_lines = dict()
real_lines = dict()
for line in lines:
    wl_cen_line = line['line_center']
    sub_region = SpectralRegion(wl_cen_line-2*u.AA, wl_cen_line+2*u.AA)
    sub_spec = extract_region(spec1D_sub, sub_region)
    real_spec = extract_region(spec1_bsmooth, sub_region)
    index = '%.1f'%wl_cen_line.value
    sub_lines[index] = sub_spec
    real_lines[index] = real_spec
#%% Plot
for key, value in sub_lines.items():
    f, ax = plt.subplots()      
    ax.step(value.spectral_axis, value.flux.value) 
#%% fitting for 5887.3
spectrum1 = sub_lines['5887.3']
g_init1 = models.Gaussian1D(amplitude=-0.0010*u.mJy, mean=5887.3*u.AA, stddev=1.*u.AA)
g_fit1 = fit_lines(spectrum1, g_init1)
y_fit1 = g_fit1(spectrum1.spectral_axis)
#%%
f, ax = plt.subplots()      
ax.step(spectrum1.spectral_axis, spectrum1.flux) 
ax.step(spectrum1.spectral_axis, y_fit1) 
#%% fitting for 5867.0
spectrum2 = sub_lines['5893.4']
g_init2 = models.Gaussian1D(amplitude=-0.0008*u.mJy, mean= 5893.4*u.AA, stddev=1.*u.AA)
g_fit2 = fit_lines(spectrum2, g_init2)
y_fit2 = g_fit2(spectrum2.spectral_axis)
f, ax = plt.subplots()      
ax.step(spectrum2.spectral_axis, spectrum2.flux) 
ax.step(spectrum2.spectral_axis, y_fit2) 
#%%

# %%
wl_D1_obs = g_fit1.mean.value
wl_D2_obs = g_fit2.mean.value

# Velocity
c = 3*10**5
v_D1 = (wl_D1_rest/wl_D1_obs - 1) * c
v_D2 = (wl_D2_rest/wl_D2_obs - 1) * c

# %% EW
spectrum1 = real_lines['5887.3']
spectrum2 = real_lines['5893.4']
ew_D2 = equivalent_width(spectrum1, continuum = 0.0025)
ew_D1 = equivalent_width(spectrum2, continuum = 0.0025)
#ew_D2 = gaussian_fwhm(spec1)
#ew_D1 = gaussian_fwhm(spec2)
f, ax = plt.subplots()      
ax.step(spec1_bsmooth.spectral_axis, spec1_bsmooth.flux, c = 'k', alpha = 0.3) 
ax.set_xlim(5880, 5900)
ax.set_ylim(0.000, 0.0035)
ax.set_xticks(np.arange(5880, 5902, 3), np.arange(5880, 5902, 3))
plt.xlabel('Wavelength[AA]')
plt.ylabel('Flux')
ax.fill_betweenx([0.00,0.0035], x1 = 5887.3-ew_D2.value/2, x2 = 5887.3+ew_D2.value/2, alpha = 0.3)
ax.step(spectrum1.spectral_axis, spectrum1.flux, label = f'Na I D2[EW = %.2f AA]'%ew_D2.value) 
ax.fill_betweenx([0.00,0.0035], x1 = 5893.4-ew_D1.value/2, x2 = 5893.4+ew_D1.value/2, alpha = 0.3)
ax.step(spectrum2.spectral_axis, spectrum2.flux, label = f'Na I D1[EW = %.2f AA]'%ew_D1.value) 
ax.legend()
#ax.step(spectrum1.spectral_axis, spectrum1.) 
# %% EW
'''
spec1 = Spectrum1D(spectral_axis = spectrum1.spectral_axis, flux = y_fit1)
spec2 = Spectrum1D(spectral_axis = spectrum2.spectral_axis, flux = y_fit2)
spec1_min_idx = np.min(np.where(np.abs(y_fit1.value) > 0.0001))
spec1_max_idx = np.max(np.where(np.abs(y_fit1.value) > 0.0001))
spec2_min_idx = np.min(np.where(np.abs(y_fit2.value) > 0.0001))
spec2_max_idx = np.max(np.where(np.abs(y_fit2.value) > 0.0001))
spec1_region = SpectralRegion(spec1.spectral_axis[spec1_min_idx], spec1.spectral_axis[spec1_max_idx])
spec2_region = SpectralRegion(spec2.spectral_axis[spec2_min_idx], spec2.spectral_axis[spec2_max_idx])
spec1 = Spectrum1D(spectral_axis = spec1D_all.spectral_axis, flux = spec1D_all.flux)
spec2 = Spectrum1D(spectral_axis = spec1D_all.spectral_axis, flux = spec1D_all.flux)
ew_D2 = equivalent_width(spec1, continuum = 1, regions = spec1_region)
ew_D1 = equivalent_width(spec2, continuum = 1, regions = spec2_region)
#ew_D2 = gaussian_fwhm(spec1)
#ew_D1 = gaussian_fwhm(spec2)
'''
#%%
spec1_bsmooth = box_smooth(spec_slice, width=30)
f, ax = plt.subplots()      
ax.step(spec1_bsmooth.spectral_axis, spec1_bsmooth.flux) 
ax.set_xlim(g_fit1.mean.value-5, g_fit1.mean.value+5)
#ax.set_xlim(g_fit2.mean.value-5, g_fit2.mean.value+5)
spec1_region = SpectralRegion(g_fit1.mean- g_fit1.stddev * 3, g_fit1.mean + g_fit1.stddev * 3)
spec2_region = SpectralRegion(g_fit2.mean- g_fit2.stddev * 3, g_fit2.mean + g_fit2.stddev * 3)
ew_D2 = equivalent_width(spec1_bsmooth, continuum = 0.00225, regions = spec1_region)
ew_D1 = equivalent_width(spec1_bsmooth, continuum = 0.00225, regions = spec2_region)
#%%
spec1_bsmooth = box_smooth(spec1D_norm, width=30)
f, ax = plt.subplots()      
ax.step(spec1_bsmooth.spectral_axis, spec1_bsmooth.flux) 
ax.set_xlim(g_fit1.mean.value-5, g_fit1.mean.value+5)
ax.set_xlim(g_fit2.mean.value-5, g_fit2.mean.value+5)
spec1_region = SpectralRegion(g_fit1.mean- g_fit1.stddev * 3, g_fit1.mean + g_fit1.stddev * 3)
spec2_region = SpectralRegion(g_fit2.mean- g_fit2.stddev * 3, g_fit2.mean + g_fit2.stddev * 3)
ew_D2 = equivalent_width(spec1_bsmooth, continuum = 1, regions = spec1_region)
ew_D1 = equivalent_width(spec1_bsmooth, continuum = 1, regions = spec2_region)

#%% plot

f, ax = plt.subplots()      
ax.set_title(r'Na ID ($D_2$)')
ax.step(spectrum1.spectral_axis, spectrum1.flux) 
#ax.step(spectrum1.spectral_axis, y_fit1)
ax.fill_betweenx([0.00,0.0028], x1 = wl_D1_obs-(ew_D2.value)/2, x2 =wl_D1_obs+(ew_D2.value)/2, alpha = 0.3)
ax.set_ylim(0.0, 0.0028)
ax.text(5885, 0.0020,r'EW = %.1f$\AA$'%ew_D2.value)
ax.text(5885, 0.0015,r'$\lambda_{cen}$ = %.1f$\AA$'%g_fit1.mean.value)

#%%
f, ax = plt.subplots()      
ax.set_title(r'Na ID ($D_1$)')
ax.step(spectrum2.spectral_axis, spectrum2.flux) 
#ax.step(spectrum1.spectral_axis, y_fit1)
ax.fill_betweenx([0.00,0.0028], x1 = wl_D2_obs-(ew_D1.value)/2, x2 =wl_D2_obs+(ew_D1.value)/2, alpha = 0.3)
ax.set_ylim(0.0, 0.0028)
ax.text(5892, 0.0020,r'EW = %.1f$\AA$'%ew_D1.value)
ax.text(5892, 0.0015,r'$\lambda_{cen}$ = %.1f$\AA$'%g_fit1.mean.value)
# %%
ebv1 = 10**(2.47*ew_D1.value - 1.76)
ebv2 = 10**(2.16*ew_D2.value - 1.91)
ew_tot = ew_D1.value + ew_D2.value
ebv_12 = 10**(1.17*ew_tot -1.85)
# %%
ebv_3 = -0.04 + 0.51 * ew_tot
# %%
ebv_3
#%%
