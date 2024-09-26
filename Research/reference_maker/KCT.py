#!/usr/bin/env python
# coding: utf-8

# In[19]:


import os, glob
import numpy as np
from File import pytable
from astropy.table import Table, QTable
import numpy as np

import os,sys,time
from astropy.io import fits
from astropy.stats import sigma_clip
import statistics
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
import pandas as pd
#%%
from astropy.coordinates import SkyCoord
from scipy import stats
from astropy.table import Table
import astroalign as aa
from astropy import units as u


# In[20]:

import sep
from Observation import cross_match
from conversion import PANSTARRS_to_JH
from Observation import to_skycoord
from Observation import to_skycoord
from Observation import cross_match
from Observation import alignIMG
from Observation import stackIMG
from conversion import SDSS_to_JH
from Observation import pytable
from File import get_fitsfilelist

from single_KCT import process_single


# In[26]:


Observatory ='KCT_STX16803'
 

# In[25]:


mode = input("Single Field/Filter or Every Field/Filter(auto) (single/auto)")


# In[27]:


option = input('choose option \n1. alignment/scaling(a/s)  \n2. Background Substraction(bkg) \n3. check seeing, depth Reference(chk) \n4. Combine aligned/scaled images(com) \n')


# In[29]:


def BKG_SUB(Field, Filter, bkgsize):
    
    inim = glob.glob(f'/data3/hhchoi1022/IMSNG/{Observatory}/selected/{Field}/{Filter}/aligned/Reference.fits')[0]
    process_single(inim)
    os.chdir(f'/data3/hhchoi1022/IMSNG/{Observatory}/selected/{Field}/{Filter}/aligned')
    
    hdr = fits.getheader(inim, ext=0)
    nim =  hdr['NCOMBINE']
    exptime = hdr['EXPTIME']
    com_exptime = int(nim*exptime)

    data = fits.getdata(inim, ext =0)
    data = data.byteswap().newbyteorder()
    bkg = sep.Background(data, bw = bkgsize, bh = bkgsize )
    substracted_data = data - bkg
    fits.writeto(f"Ref-{Observatory}-{Field}-{Filter}-{com_exptime}.com.fits", data, hdr, overwrite = True)
    fits.writeto(f"Ref-subtracted-{Observatory}-{Field}-{Filter}-{com_exptime}.com.fits", substracted_data, hdr, overwrite = True)


# In[ ]:
def SCALE(image,tgtim):
    tgtim_data = fits.getdata(tgtim, ext=0)
    tgtim_data = tgtim.byteswap().newbyteorder()
    tgtim_hdr = fits.getheader(tgtim)
    image_data = fits.getdata(image,ext=0)
    image_data = image_data.byteswap().newbyteorder()
    image_hdr = fits.getheader(image)
    bkg_offset = image_hdr['SKYVAL'] - tgtim_hdr['SKYVAL']
    scaled_data = image_data-bkg_offset
    return scaled_data, image_hdr

def ALIGN(image, refimage):
    reference_data = fits.getdata(refimage, ext=0)
    reference_data = reference_data.byteswap().newbyteorder()
    reference_hdr = fits.getheader(refimage)
    image_data, image_hdr = SCALE(image,refimage)
    aligned_data, footprint = aa.register(image_data, reference_data)
    return aligned_data, reference_hdr
    

def ALIGN_SCALE(Field, Filter):
    path = f'/data3/hhchoi1022/IMSNG/{Observatory}/selected/{Field}/{Filter}'
    os.chdir(path)

    imkey = f'{path}/*.fits'
    imlist = sorted(glob.glob(f'{imkey}'))
    outbl = Table()
    outbl['Image'] = imlist
    outbl['Seeing'] = 0.0
    outbl['Depth'] = 0.0
    outbl['ZP'] = 0.0

    def gethdrinfo(inim, key):
        hdr = fits.getheader(inim)
        if key in hdr.keys():
            value = hdr[key]
        else:
            value = None
        return value
    # imlist = removed image list
    
    for i, inim in enumerate(imlist):
        outbl['Seeing'][i] = gethdrinfo(inim, 'HSEEING')
        outbl['Depth'][i] = gethdrinfo(inim, 'HDEPTH5')
        outbl['ZP'][i] = gethdrinfo(inim,'HZP_AP')
    os.makedirs('aligned', exist_ok = True)
    os.chdir('aligned')

    os.getcwd()
    reference_image = outbl['Image'][np.argmax(outbl['ZP'])]
    print(f"reference image path = {os.path.basename(reference_image)}")

    for image in imlist:
        aligned_data, reference_hdr = ALIGN(image,reference_image)


    # images -> scaled image aligned for image combine 
    for i, inim in enumerate(imlist):
        try:
            start = time.perf_counter()
            hdr = fits.getheader(inim, ext=0)
            data = fits.getdata(inim, ext =0)
            data = data.byteswap().newbyteorder()
            #bkg = sep.Background(data, bw = bkgsize, bh = bkgsize )
            substracted_data = data# - bkg
            scaled_substracted_data = substracted_data*100**(-DELMAG[i]/5)
            aligned_data, footprint = aa.register(scaled_substracted_data, reference_data)
            fits.writeto(f"{os.path.basename(inim)}", aligned_data, reference_hdr)
            end = time.perf_counter()-start
            print(f'{i}/{len(imlist)} image processed \nExpected time remaining : ', int(end*(len(imlist)-i)),'sec')
        except :
            print(f'{inim} : Alignment Failed')




#%%
def COMBINE(target_, filter_):
    import os
    
    os.chdir('/data2/iraf')
    
    from pyraf import iraf
    import glob
    from astropy.io import fits
    from astropy.io.fits import getheader
    #from single_KCT import process_single
    
    targetlist = os.listdir('/data3/hhchoi1022/IMSNG/KCT_STX16803/selected/')
    targetlist = targetlist[19:]
    filterlist = ['g','r','i']
    refdir = '/data3/hhchoi1022/IMSNG/KCT_STX16803/Reference'

    imdir = f'/data3/hhchoi1022/IMSNG/KCT_STX16803/selected/{target_}/{filter_}/aligned'
    imkey = f'{imdir}/C*.fits'
    imagelist = glob.glob(imkey)
    n_im = len(imagelist)
    if not n_im == 0:

        hdr = fits.getheader(imagelist[0])
        exptime = hdr['EXPTIME']
        tot_exptime = int(exptime*n_im)
        
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
        


# In[ ]:


r_depth = []
r_seeing =[]
c_depth = []
c_seeing =[]

Observatory = 'KCT_STX16803'
PIXSCALE = 0.724
Filterset = ['g','r','i'] ##########FILTER##########

if mode == 'single':
    try :
        Field = input("Input Field name")
        Filterlist = ['g','r','i']
        #input("Input Filter name")
        Field_Catalog =  pytable(f'/data2/hhchoi1022/IMSNG/alltarget.dat') 
        pixsize=Field_Catalog[Field_Catalog['obj']==Field]['maxaxis']*60/PIXSCALE   
        
        if pixsize <=64:
            bkgsize = 64
        elif pixsize > 64 and pixsize <=128:
            bkgsize = 128
        elif pixsize >128 and pixsize <=400:
            bkgsize = 256  
        elif pixsize >400:
            bkgsize = 512  
        for Filter in Filterlist:    
            if option == 'a/s':
                ALIGN_SCALE(Field, Filter)
            elif option == 'bkg':  
                BKG_SUB(Field, Filter,bkgsize)
            elif option == 'com':  
                COMBINE(Field, Filter)
            elif option == 'auto':
                ALIGN_SCALE(Field, Filter)
                COMBINE(Field, Filter)
                BKG_SUB(Field, Filter,bkgsize)
        
        os.system(f'cp {path}/aligned/Ref-*.fits /data3/hhchoi1022/IMSNG/{Observatory}/Reference/')

    except:
        print("__________WRONG value____________")
elif mode == 'auto' or mode == '':
    Field_Catalog =  pytable(f'/data2/hhchoi1022/IMSNG/alltarget.dat')
    fieldlist = os.listdir('/data3/hhchoi1022/IMSNG/KCT_STX16803/selected/')
    Catalog = []
    for Field in fieldlist:
        pixsize=Field_Catalog[Field_Catalog['obj']==Field]['maxaxis']*60/PIXSCALE
        Catalog.append([Field,int(pixsize)])
    
    for i, (Field, Angular_size) in enumerate(Catalog):
        if Angular_size <=64:
            bkgsize = 64
        elif Angular_size > 64 and Angular_size <=128:
            bkgsize = 128
        elif Angular_size >128 and Angular_size <=400:
            bkgsize = 256  
        elif Angular_size >400:
            bkgsize = 512  
        for Filter in Filterset:
            path = f'/data3/hhchoi1022/IMSNG/{Observatory}/selected/{Field}/{Filter}'
            try:
                if os.path.isdir(path):
                    if not os.path.isdir(f'{path}/aligned'):
                        print(f'aligning {Field}-{Filter}')
                        #ALIGN_SCALE(Field,Filter)
                    if not os.path.isfile(f'{path}/aligned/Reference.fits'):
                        print(f'Combining {Field}-{Filter}')
                        #COMBINE(Field,Filter)
                        BKG_SUB(Field, Filter,bkgsize)

                        print(90*'-')
                        os.system(f'cp {path}/aligned/Ref-*.fits /data3/hhchoi1022/IMSNG/{Observatory}/Reference/')
                    if os.path.isfile(f'{path}/aligned/Reference.fits'):
                        BKG_SUB(Field, Filter,bkgsize)

                        print(90*'-')
                        os.system(f'cp {path}/aligned/Ref-*.fits /data3/hhchoi1022/IMSNG/{Observatory}/Reference/')

            except:
                print(f'{Field}-{Filter} cannot be calculated')

                
#                 if os.path.isdir(f'{path}/aligned'):
#                     if option == 'bkg':
#                         if os.path.isfile(f'{path}/aligned/Reference.fits'):
#                             Refbkg_subs(Field, Filter,bkgsize)
#                             os.system(f'cp {path}/aligned/Ref-*.fits /data3/hhchoi1022/IMSNG/{Observatory}/Reference/')

#                     if option == 'chk':
#                         if os.path.isfile(f'{path}/aligned/Reference.fits'):
#                             Ref_outbl, exp_depth, exp_seeing = check_IMG_from_hdr(Field,Filter)
#                             Cal_outbl = check_depth_seeing(Field, Filter,bkgsize)
#                             plt.plot(Ref_outbl['Seeing'],Ref_outbl['Depth'], marker = 'o', mec=  'k', mfc = 'none', ls = 'none', alpha = 0.25)
#                             plt.plot(Cal_outbl['Seeing'],Cal_outbl['Depth'], marker = 'o', mec=  'red', mfc = 'none', ls = 'none', alpha = 0.55)
#                             plt.plot(exp_seeing, exp_depth, marker = '+', color = 'blue', alpha = 0.5)
#                             plt.xlabel('Seeing [arcsec]', fontsize=20)
#                             plt.ylabel('Depth [AB]', fontsize=20)
#                             plt.grid('both', ls='--', c='grey', alpha=0.5)
#                             plt.legend()
#                             plt.tight_layout()
#                             plt.minorticks_on()
#                             plt.show()
#                             r_depth.append(exp_depth)
#                             r_seeing.append(exp_seeing)
#                             c_depth.append(Cal_outbl['Depth'])
#                             c_seeing.append(Cal_outbl['Seeing'])
#                     if option == 'com':
#                         if not os.path.isfile(f'{path}/aligned/Reference.fits'):
#                             print(f'{Field}/{Filter} processing...')
#                             start=  time.perf_counter()
#                             IMCOMBINE(Field,Filter)
#                             print(61*(time.perf_counter()-start),'sec remaining')
#                         else:
#                             print(f'{Field}-{Filter}-Reference already exists')
#                 else:
#                     if option == 'a/s':
#                         align_scaling(Field,Filter)

# In[ ]:


