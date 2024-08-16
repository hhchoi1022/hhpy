#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 19:52:40 2021

@author: hhchoi10222
"""

import os

os.chdir('/data2/iraf')

from pyraf import iraf
import glob
from astropy.io import fits
from astropy.io.fits import getheader
#from single_KCT import process_single
#%%
targetlist = os.listdir('/data3/hhchoi1022/IMSNG/KCT_STX16803/selected/')
targetlist = targetlist[19:]
filterlist = ['g','r','i']
refdir = '/data3/hhchoi1022/IMSNG/KCT_STX16803/Reference'
for target_ in targetlist:
    for filter_ in filterlist:
        imdir = f'/data3/hhchoi1022/IMSNG/KCT_STX16803/selected/{target_}/{filter_}/aligned'
        imkey = f'{imdir}/C*.fits'
        imagelist = glob.glob(imkey)
        n_im = len(imagelist)
        if n_im == 0:
            continue
        hdr = fits.getheader(imagelist[0])
        exptime = hdr['EXPTIME']
        tot_exptime = int(exptime*n_im)
        imdir = '/data2/temporary'
        os.chdir(imdir)
        iraf.chdir(imdir)
        os.system('ls C*.fits > image.list')
        if n_im > 30:
            iraf.imcombine('@image.list', 'Reference.fits', combine = 'median', reject = 'minmax', zero = 'mode',nlow = 2, nhigh= 4)
        elif n_im > 15:
            iraf.imcombine('@image.list', 'Reference.fits', combine = 'median', reject = 'minmax', zero = 'mode',nlow = 1, nhigh= 3)
        else:
            iraf.imcombine('@image.list', 'Reference.fits', combine = 'median', reject = 'minmax', zero = 'mode',nlow = 1, nhigh= 1)

        os.system('rm image.list')
        
        os.system(f'cp Reference.fits {refdir}')
        os.chdir(refdir)
        os.system(f'mv Reference.fits Ref-KCT_STX16803-{target_}-{filter_}-{tot_exptime}.com.fits')
        
        

