#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 14:22:28 2021

@author: hhchoi10222
"""
from astropy.io import fits
import sep
import numpy as np
import matplotlib.pyplot as plt
from astropy.io.fits import getdata
import os




import glob
#%%

def bkgqual(imagelist):
    bkg_std = []
    bkg_val = []
    for i,image in enumerate(imagelist):
        data = fits.getdata(image)
        data = data.byteswap().newbyteorder()
        bkg = sep.Background(data)
        bkg_std.append(bkg.globalrms)
        bkg_val.append(bkg.globalback)
        
        plt.figure()
        plt.title('bkg map')
        plt.imshow(bkg)    
    index = np.where(bkg_std>np.median(bkg_std)+1.5*np.std(bkg_std))
    badimlist = np.array(imagelist)[index]
    return badimlist


fieldlist = os.listdir('/data3/hhchoi1022/IMSNG/KCT_STX16803/selected/')
filterlist = ['g','r','i']
info = []
for field_ in fieldlist:
    for filter_ in filterlist:
        
        print(60*'-')
        print(field_, filter_)
        imkey = f'/data3/hhchoi1022/IMSNG/KCT_STX16803/selected/{field_}/{filter_}/*.fits'
        imagelist = sorted(glob.glob(imkey))
        command = 'rm'
        if len(imagelist) > 2:
            b_count = len(imagelist)
            badimlist= bkgqual(imagelist)
            rm_count = len(badimlist)
            info.append([field_,b_count,rm_count])
            
            print(rm_count, b_count)
            if rm_count/b_count < 0.3:
                
                for imname in badimlist:
                    command = command+' '+imname
                print(command)
                #os.system(command)
            else:
                if b_count > 30:
                    for imname in badimlist:
                        command = command+' '+imname
                    print(command)
        if not command =='rm':
            print(command)
            os.system(command)
            
                    

            
    