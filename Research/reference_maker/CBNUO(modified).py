#!/usr/bin/env python
# coding: utf-8

# In[10]:


import os, glob
import numpy as np
from File import pytable
from astropy.table import Table, QTable
import astroalign as aa
import numpy as np
import sep
import time
from astropy.io import fits
from astropy.stats import sigma_clip
import statistics
import matplotlib.pyplot as plt


# In[11]:


from Observation import cross_match
from conversion import PANSTARRS_to_JH
from Observation import to_skycoord


# In[12]:


from File import get_fitsfilelist


# In[13]:


from astropy.nddata import CCDData


# In[14]:


from astropy import units as u


# In[15]:


from ccdproc import Combiner


# In[16]:


Observatory = 'CBNUO'


# In[17]:


mode = input("Single Field/Filter or Every Field/Filter(auto) (single/auto)")


# In[18]:


option = input('choose option \n1. alignment/scaling(a/s)  \n2. Background Substraction(bkg) \n3. check seeing, depth Reference(chk) \n4. Combine aligned/scaled images(com) \n')


# In[10]:


from Single_KCT import process_single


# In[13]:


def Refbkg_subs(Field, Filter, bkgsize):
    
    inim = glob.glob(f'/data3/hhchoi1022/IMSNG/{Observatory}/Selected/{Field}/{Filter}/aligned/Reference.fits')[0]
    process_single(inim)
    os.chdir(f'/data3/hhchoi1022/IMSNG/{Observatory}/Selected/{Field}/{Filter}/aligned')
    
    hdr = fits.getheader(inim, ext=0)
    nim =  hdr['NCOMBINE']
    exptime = hdr['EXPTIME']
    com_exptime = int(nim*exptime)

    data = fits.getdata(inim, ext =0)
    data = data.byteswap().newbyteorder()
    bkg = sep.Background(data, bw = bkgsize, bh = bkgsize )
    substracted_data = data - bkg
    fits.writeto(f"Ref-CBNUO-{Field}-{Filter}-{com_exptime}.com.fits", data, hdr, overwrite = True)
    fits.writeto(f"Ref-substracted-CBNUO-{Field}-{Filter}-{com_exptime}.com.fits", substracted_data, hdr, overwrite = True)


# In[ ]:





# In[14]:


def align_scaling(Field, Filter):
    path = f'/data3/hhchoi1022/IMSNG/{Observatory}/Selected/{Field}/{Filter}'
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
    reference_data = fits.getdata(reference_image, ext=0)
    reference_data = reference_data.byteswap().newbyteorder()
    reference_hdr = fits.getheader(reference_image)
    print(f"reference image path = {reference_image}")
    DELMAG = (outbl['ZP']-np.max(outbl['ZP'])).round(4)


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
            print(f'{i} image processed \nExpected time remaining : ', int(end*(len(imlist)-i)),'sec')
        except :
            print(f'{inim} : Alignment Failed')


# In[33]:


def IMCOMBINE(Field, Filter):
    path = f'/data3/hhchoi1022/IMSNG/{Observatory}/Selected/{Field}/{Filter}/aligned'
    os.chdir(path)
    imkey = glob.glob(f'{path}/C*.fits')
    
    DATA = []
    for image in imkey:
        hsulist = fits.open(image)
        data = hsulist[0].data
        bkg = np.median(data)
        data = data-bkg
        ccddata = CCDData(data, unit = u.adu)
        DATA.append(ccddata)
    
    #header 
    hdr= fits.getheader(imkey[0])
    hdr.append('COMATHR', end = True)
    hdr['COMATHR'] = ('Hyeonho Choi')
    for i, image in enumerate(imkey):
        hdr.append(f'IMCOM{i}', end = True)
        hdr[f'IMCOM{i}'] = (os.path.basename(image))
    hdr.append('NCOMBINE', end =True)
    hdr['NCOMBINE'] = len(imkey)
    
    #Combine
    combiner = Combiner(DATA)
    combined_median = combiner.median_combine()
    
    fits.writeto('Reference.fits',combined_median, hdr)


# In[34]:


r_depth = []
r_seeing =[]
c_depth = []
c_seeing =[]
if mode == 'single':
    try :
        if option == 'a/s':
            Field = input("Input Field name")
            Filter = input("Input Filter name")
            align_scaling(Field, Filter)
        elif option == 'bkg':
            Field = input("Input Field name")
            Filter = input("Input Filter name")
            bkgsize = int(input('bkgsize'))   
            Refbkg_subs(Field, Filter,bkgsize)
        elif option == 'com':
            Field = input("Input Field name")
            Filter = input("Input Filter name")   
            IMCOMBINE(Field, Filter)
    except:
        print("__________WRONG value____________")
elif mode == 'auto' or mode == '':
    Field_Catalog = input('Input Field Catalog path(*/Field_Catalog) or default(default)')
    if Field_Catalog == 'default' or Field_Catalog == '':
        Field_Catalog = f'/data3/hhchoi1022/IMSNG/{Observatory}/Field_Catalog'
    Catalog = pytable(Field_Catalog)
    Filterset = ['B','V','R','I']
    for i, (Field, RA, Dec, Angular_size) in enumerate(Catalog):
        if Angular_size <= 64:
            bkgsize = 64
        elif Angular_size > 64 and Angular_size <=192:
            bkgsize = 128
        elif Angular_size >192 :
            bkgsize = 256      
        for Filter in Filterset:
            path = f'/data3/hhchoi1022/IMSNG/{Observatory}/Selected/{Field}/{Filter}'
            if os.path.isdir(path):
                if os.path.isdir(f'{path}/aligned'):
                    if option == 'bkg':
                        if os.path.isfile(f'{path}/aligned/Reference.fits'):
                            Refbkg_subs(Field, Filter,bkgsize)
                            os.system(f'cp {path}/aligned/Ref-*.fits /data3/hhchoi1022/IMSNG/{Observatory}/Reference/')

                    if option == 'chk':
                        if os.path.isfile(f'{path}/aligned/Reference.fits'):
                            Ref_outbl, exp_depth, exp_seeing = check_IMG_from_hdr(Field,Filter)
                            Cal_outbl = check_depth_seeing(Field, Filter,bkgsize)
                            plt.plot(Ref_outbl['Seeing'],Ref_outbl['Depth'], marker = 'o', mec=  'k', mfc = 'none', ls = 'none', alpha = 0.25)
                            plt.plot(Cal_outbl['Seeing'],Cal_outbl['Depth'], marker = 'o', mec=  'red', mfc = 'none', ls = 'none', alpha = 0.55)
                            plt.plot(exp_seeing, exp_depth, marker = '+', color = 'blue', alpha = 0.5)
                            plt.xlabel('Seeing [arcsec]', fontsize=20)
                            plt.ylabel('Depth [AB]', fontsize=20)
                            plt.grid('both', ls='--', c='grey', alpha=0.5)
                            plt.legend()
                            plt.tight_layout()
                            plt.minorticks_on()
                            plt.show()
                            r_depth.append(exp_depth)
                            r_seeing.append(exp_seeing)
                            c_depth.append(Cal_outbl['Depth'])
                            c_seeing.append(Cal_outbl['Seeing'])
                    if option == 'com':
                        if not os.path.isfile(f'{path}/aligned/Reference.fits'):
                            print(f'{Field}/{Filter} processing...')
                            start=  time.perf_counter()
                            IMCOMBINE(Field,Filter
                                     )
                            print(61*(time.perf_counter()-start),'sec remaining')
                        else:
                            print(f'{Field}-{Filter}-Reference already exists')
                else:
                    if option == 'a/s':
                        align_scaling(Field,Filter)


# In[ ]:




