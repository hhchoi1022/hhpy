#%%
import astropy.units as u
from specutils import Spectrum1D
from specutils.manipulation import (box_smooth, gaussian_smooth, trapezoid_smooth)
#import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
# %%
'''
tar = ascii.read('/home/hhchoi1022/Desktop/2022-11-04/pNGC2392_w_spec.csv')
tar_wave = tar['wavelength']
tar_cal = tar['flux']
tar_std = tar['fluxerr']
fig,ax = plt.subplots(1,1,figsize=(10,6))
ax.plot(tar_wave,tar_cal,c='k')
ax.set_xlim(4000,7500)
ax.set_title('After flux calibration')
ax.set_ylabel(r'erg $s^{-1}cm^{-2}\AA^{-1}$ ')
ax.set_xlabel(r'Dispersion axis $[\AA]$',fontsize=20) 
#fig.savefig(SUBPATH/('p'+path.stem+'_w_spec.png'))
# %%
def spectrum(spec_tbl):
    spec = Spectrum1D(flux = spec_tbl['flux']*u.Unit('erg cm-2 s-1 AA-1') , spectral_axis = spec_tbl['wavelength']*u.AA)
    return spec
# %%


spec_smooth = box_smooth(spectrum(tar), width  =40)
#%%
fig,ax = plt.subplots(1,1,figsize=(10,6))
ax.plot(tar_wave,tar_cal,c='k')
ax.plot(spec_smooth.spectral_axis , spec_smooth.flux,c='r', linestyle= '--')
ax.set_xlim(4000,7500)
ax.set_title('After flux calibration')
ax.set_ylabel(r'erg $s^{-1}cm^{-2}\AA^{-1}$ ')
ax.set_xlabel(r'Dispersion axis $[\AA]$',fontsize=20) 
'''
# %%
import pylab as plt
import numpy as np
from scipy.interpolate import splrep,splev
import sys
import os

def onclick(event):
    # when none of the toolbar buttons is activated and the user clicks in the
    # plot somewhere, compute the median value of the spectrum in a 10angstrom
    # window around the x-coordinate of the clicked point. The y coordinate
    # of the clicked point is not important. Make sure the continuum points
    # `feel` it when it gets clicked, set the `feel-radius` (picker) to 5 points
    toolbar = plt.get_current_fig_manager().toolbar
    if event.button==1 and toolbar.mode=='':
        window = ((event.xdata-5)<=wave) & (wave<=(event.xdata+5))
        y = np.median(flux[window])
        plt.plot(event.xdata,y,'rs',ms=10,picker=5,label='cont_pnt')

    plt.draw()

def onpick(event):
    # when the user clicks right on a continuum point, remove it
    if event.mouseevent.button==3:
        if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
            event.artist.remove()

def ontype(event):
    # when the user hits enter:
    # 1. Cycle through the artists in the current axes. If it is a continuum
    #    point, remember its coordinates. If it is the fitted continuum from the
    #    previous step, remove it
    # 2. sort the continuum-point-array according to the x-values
    # 3. fit a spline and evaluate it in the wavelength points
    # 4. plot the continuum
    if event.key=='enter':
        cont_pnt_coord = []
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='cont_pnt':
                cont_pnt_coord.append(artist.get_data())
            elif hasattr(artist,'get_label') and artist.get_label()=='continuum':
                artist.remove()
        cont_pnt_coord = np.array(cont_pnt_coord)[...,0]
        sort_array = np.argsort(cont_pnt_coord[:,0])
        x,y = cont_pnt_coord[sort_array].T
        spline = splrep(x,y,k=3)
        continuum = splev(wave,spline)
        cont_tbl = Table()
        cont_tbl['wave'] = wave
        cont_tbl['flux'] = continuum
        cont_tbl.write('NGC2392_cont.dat', format = 'ascii', overwrite= True)
        plt.plot(wave,continuum,'r-',lw=2,label='continuum')
        

    # when the user hits 'n' and a spline-continuum is fitted, normalise the
    # spectrum
    elif event.key=='n':
        continuum = None
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='continuum':
                continuum = artist.get_data()[1]
                break
        if continuum is not None:
            plt.cla()
            plt.plot(wave,flux/continuum,'k-',label='normalised')

    # when the user hits 'r': clear the axes and plot the original spectrum
    elif event.key=='r':
        plt.cla()
        plt.plot(wave,flux,'k-')

    # when the user hits 'w': if the normalised spectrum exists, write it to a
    # file.
    elif event.key=='w':
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='normalised':
                data = np.array(artist.get_data())
                np.savetxt(os.path.splitext(filename)[0]+'.nspec',data.T)
                print('Saved to file')
                break
    plt.draw()
#%%

if __name__ == "__main__":
    # Get the filename of the spectrum from the command line, and plot it
    filename = sys.argv[1]
    data = ascii.read(filename)
    wave, flux = data['wavelength'], data['flux']
    spectrum, = plt.plot(wave,flux,'k-',label='spectrum')
    plt.title(filename)

    # Connect the different functions to the different events
    plt.gcf().canvas.mpl_connect('key_press_event',ontype)
    plt.gcf().canvas.mpl_connect('button_press_event',onclick)
    plt.gcf().canvas.mpl_connect('pick_event',onpick)
    plt.show() # show the window
# %%
filename = './pNGC2392_w_spec.csv'
data = ascii.read(filename)
wave, flux = data['wavelength'], data['flux']
spectrum, = plt.plot(wave,flux,'k-',label='spectrum')
plt.title(filename)
plt.figure(dpi = 600)
plt.plot(wave, flux)

# Connect the different functions to the different events
plt.gcf().canvas.mpl_connect('key_press_event',ontype)
plt.gcf().canvas.mpl_connect('button_press_event',onclick)
plt.gcf().canvas.mpl_connect('pick_event',onpick)


# %%
#%%
filename = './pNGC2392_w_spec.csv'
data = ascii.read(filename)
wave, flux = data['wavelength'], data['flux']
spectrum, = plt.plot(wave,flux,'k-',label='spectrum')
# %%
continuum_filename = './NGC2392_cont.dat'
data_cont = ascii.read(continuum_filename)
wave_cont, flux_cont = data_cont['wave'], data_cont['flux']
plt.plot(wave_cont,flux - flux_cont,'k-',label='spectrum')
# %%

# %%
wave