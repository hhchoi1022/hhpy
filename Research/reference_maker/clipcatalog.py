#!/usr/bin/env python


# In[20]:


import os,glob
import csv
from astropy.io import ascii
import pandas as pd


# In[2]:


from Observation import hms2dec, dms2dec


# In[9]:


from File import pytable


# In[12]:


targetlist = '/data2/hhchoi1022/IMSNG/alltarget.dat'
target = pytable(targetlist)


# In[27]:


def clipping_catalog(catalog, fov):
    table = pd.read_csv(catalog)
    targetname = os.path.basename(catalog).split('_')[0]
    target_info = target[target['obj']==targetname]
    RA = target_info['ra'][0]
    Dec = target_info['dec'][0]
    center = [hms2dec(RA),dms2dec(Dec)]
    FOV = fov/60
    
    ## variables
    os.chdir('/data3/hhchoi1022/IMSNG/KCT_STX16803/SkyCatalog/Skymapper/')
    clipped_table = table[(table['raj2000']<center[0]+FOV)&(table['raj2000']>center[0]-FOV)&(table['dej2000']<center[1]+FOV)&(table['dej2000']>center[1]-FOV)]
    ####
    try:
        clipped_table.to_csv(f'{targetname}.csv',
                    sep = ',',
                    na_rep ='',
                    index =False)
        print(f'{targetname} Succelfully downloaded')
    except:
        print(f'{targetname} went wrong')


# In[16]:


imkey = '/data2/hhchoi1022/IMSNG/Skymapper/*.csv'


# In[29]:


import numpy as np


# In[23]:


catlist = glob.glob(imkey)


# In[28]:


for catalog in catlist:
    try:
        clipping_catalog(catalog,60)
    except:
        print('clipping failed')


`# In[ ]:




