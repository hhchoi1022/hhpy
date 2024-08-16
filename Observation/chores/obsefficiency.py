#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 10:22:18 2022

@author: hhchoi1022
"""

from astropy.io import fits
import os, glob
from astropy.time import Time
from datetime import datetime
from matplotlib import dates
import numpy as np
import matplotlib.pyplot as plt
date = '20220404'

datelist = sorted(os.listdir('/data2/temp/'))

timelist1 = []
timelist2 = []
timelist3= []
yaxis = []
ratio1 = []
for date in datelist:
    dirname = f'/data2/temp/{date}'
    imkey = f'{dirname}/*.fts'
    imlist = sorted(glob.glob(imkey))
    timelist = []
    for image in imlist:
        hdr = fits.getheader(image)
        time = Time(hdr['DATE-OBS']).iso
        datet = datetime.strptime(time, '%Y-%m-%d %H:%M:%S.%f')
        datey = int(datetime.strptime(date, '%Y%m%d').day)
        yaxis.append(datey)
        timelist.append(datet)
    timelist.sort()
    timelist1 = [(timelist[i]-timelist[0]).seconds for i in range(len(timelist))]
    timelist2 = np.array(timelist1)/3600
    timelist3.append(timelist2)
    
    sumtime = 0
    for i in range(len(timelist2)):
        if (timelist2[i]-timelist2[i-1]) > 0.10:
            sumtime = sumtime+(timelist2[i]-timelist2[i-1])
    tottime = np.abs(timelist2[0]-timelist2[-1])
    ratio = round((1-sumtime/tottime)*100,1)
    ratio1.append(ratio)

timelist4 =[]
for timelist in timelist3:
    for time in timelist:
        timelist4.append(time)
#%%
plt.figure(dpi = 300)
plt.xlim(-1,10)
plt.xlabel('Hours')
plt.ylabel('Days[4/?]')
plt.plot(timelist4, yaxis, marker= '.',linewidth = 0, c = 'k', label = ratio1)
#plt.legend()
plt.xticks(rotation=45)

