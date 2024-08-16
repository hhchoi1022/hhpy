#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 09:58:12 2022

@author: hhchoi1022
"""
from astropy.io import fits
import numpy as np
import glob
import os
import sys
import os, glob
from astropy.io import fits


imkey = '/data2/SN2021aefx_220627/RASA36/r/NOINFO/*.fits'
imlist = sorted(glob.glob(imkey))
bkglist = []


def divide_RASA(imlist):
    for i, image in enumerate(imlist):
        printProgress(i,len(imlist),prefix = 'Unsorted',suffix = 'Sorted')
        imdir = os.path.dirname(image)
        hdr = fits.getheader(image)
        if 'CAMMODE' in hdr.keys():
            if hdr['CAMMODE'] == 'HIGH':
                os.makedirs(f'{imdir}/HIGH', exist_ok = True)
                os.system(f'cp {image} {imdir}/HIGH')
            else:
                os.makedirs(f'{imdir}/MERGE', exist_ok = True)
                os.system(f'cp {image} {imdir}/MERGE')
        else:
            os.makedirs(f'{imdir}/NOINFO', exist_ok = True)
            os.system(f'cp {image} {imdir}/NOINFO')

def divide_value(imlist):
    for i, image in enumerate(imlist):
        printProgress(i,len(imlist),prefix = 'Unsorted',suffix = 'Sorted')
        imdir = os.path.dirname(image)
        data = fits.getdata(image)
        zero = np.mean(np.sort((data[~np.isnan(data)].flatten()))[-10:-5])
        bkglist.append(zero)
        if zero > 80000:
            os.makedirs(f'{imdir}/MERGE', exist_ok = True)
            os.system(f'mv {image} {os.path.dirname(image)}/MERGE')
            print(f'{image} is sorted as MERGE')
        elif zero < 80000:
            os.makedirs(f'{imdir}/HIGH', exist_ok = True)
            os.system(f'mv {image} {os.path.dirname(image)}/HIGH')
            print(f'{image} is sorted as HIGH')

def printProgress (iteration, total, prefix = '', suffix = '', decimals = 1, barLength = 100): 
    formatStr = "{0:." + str(decimals) + "f}" 
    percent = formatStr.format(100 * (iteration / float(total))) 
    filledLength = int(round(barLength * iteration / float(total))) 
    bar = '#' * filledLength + '-' * (barLength - filledLength) 
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percent, '%', suffix)), 
    if iteration == total: 
        sys.stdout.write('\n') 
    sys.stdout.flush() 
#%%
imlist
divide_value(imlist)







    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    