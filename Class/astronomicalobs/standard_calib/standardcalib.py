#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 15:38:57 2022

@author: hhchoi1022
"""

# Photometric calibration 

from astropy.table import Table
import numpy as np
tbl = Table()

objectlist = ['a1', 'a2', 'a3', 'a4', 'a5'] *3 + ['b1', 'b2', 'b3', 'b4', 'b5'] *2 + ['s1','s2','q1','g1']
airmasslist = [1]*5+[1.3]*5+[2]*5+[1.5]*5+[2.5]*5+[1.2,1.2,1.4,2.1]
expVlist = [5]*5+[6]*5+[7]*5+[10]*5+[5]*5+[15,15,200,600]
expBlist = [10]*5+[12]*5+[15]*5+[20]*5+[10]*5+[30,30,400,1200]
Vcountlist = [12687301, 5534312, 2280320, 990631, 337942, 14711983, 6404748, 2638025, 1138554, 388549, 15845201, 6885820, 2796503, 1206461, 410667, 16611138, 7207806, 2238195, 1400365, 332610, 7395540, 3208276, 979923, 610134, 142919, 8695591, 1998330, 180181, 1556]
Bcountlist = [22438244, 7521077, 1835230, 560093, 124145, 25343232, 8489650, 2060068, 626563, 138593, 25669738, 8636456, 2042511, 630063, 130178, 25755872, 9408614, 1881269, 584810, 74996, 10503907, 3887342, 740349, 236381, 28557, 14211249, 754762, 112201, 835]

tbl['object'] = objectlist
tbl['airmass'] = airmasslist
tbl['expV'] = expVlist
tbl['expB'] = expBlist
tbl['Vcountlist'] = Vcountlist
tbl['Bcountlist'] = Bcountlist
tbl['Vcountlist_1s'] = tbl['Vcountlist']/tbl['expV']
tbl['Bcountlist_1s'] = tbl['Bcountlist']/tbl['expB']
tbl['Berror'] = np.sqrt(tbl['Vcountlist'])
tbl['Verror'] = np.sqrt(tbl['Bcountlist'])
tbl['Berror_1s'] = np.sqrt(tbl['Bcountlist_1s'])
tbl['Verror_1s'] = np.sqrt(tbl['Vcountlist_1s'])
tbl.write('tbl.dat', format = 'ascii')
tbl_ref = Table()

objectlist = ['a1', 'a2', 'a3', 'a4', 'a5'] *1 + ['b1', 'b2', 'b3', 'b4', 'b5'] *1
vlist = [9.1, 10.0, 11.0, 11.9, 13.1, 9.5, 10.4, 11.7, 12.2, 13.8]
BVlist = [-0.3, 0, 0.6, 1.0, 1.5, -0.2, 0, 0.5, 1.3, 2.0]

tbl_ref['object'] = objectlist
tbl_ref['B-V'] = BVlist
tbl_ref['VMAG'] = vlist
tbl_ref['BMAG'] = tbl_ref['B-V'] + tbl_ref['VMAG']

#tbl = join(tbl_ref, tbl, keys = 'object')

#%%

import numpy as np
import matplotlib.pyplot as plt

tbl['Vmag'] = -2.5*np.log10(tbl['Vcountlist_1s'])+25
tbl['Bmag'] = -2.5*np.log10(tbl['Bcountlist_1s'])+25
tbl['Vmagerr'] = 2.5*tbl['Verror']/2.303/tbl['Vcountlist']
tbl['Bmagerr'] = 2.5*tbl['Berror']/2.303/tbl['Bcountlist']
tbl['BVmag'] = tbl['Bmag'] - tbl['Vmag']

#%% k'
from astropy.table import vstack
group_obj = tbl.group_by('object')
a2group = group_obj.groups[1]
tbl1 = tbl[(tbl['object'] == 'a2')]
tbl1_ref = tbl_ref[tbl_ref['object'] == 'a2'] 
a2group['VMAG'] = tbl1_ref['VMAG']
a2group['BMAG'] = tbl1_ref['B-V']+tbl1_ref['VMAG']
a2group['B-V'] = tbl1_ref['B-V']
b2group = group_obj.groups[6]
tbl2 = tbl[(tbl['object'] == 'b2')]
tbl2_ref = tbl_ref[tbl_ref['object'] == 'b2']
b2group['VMAG'] = tbl2_ref['VMAG']
b2group['BMAG'] = tbl2_ref['VMAG']
b2group['B-V'] = tbl1_ref['B-V']+tbl2_ref['B-V']
ab2group = vstack([a2group, b2group])

airlist = ab2group['airmass']
Vmaglist = ab2group['Vmag']-ab2group['VMAG']
Bmaglist = ab2group['Bmag']-ab2group['BMAG']

from scipy.optimize import curve_fit

def func(x, a, b):
    y = a*x + b
    return y

plt.figure(figsize = (6,4), dpi = 500)
x = np.linspace(np.min(airlist),np.max(airlist))
popt, pcov = curve_fit(func, airlist, Vmaglist, sigma = ab2group['Vmagerr'])
plt.errorbar(x = a2group['airmass'],y = a2group['Vmag']-a2group['VMAG'], yerr = a2group['Vmagerr'], fmt = 'D', elinewidth  = 2, c = 'k')
plt.errorbar(x = b2group['airmass'],y = b2group['Vmag']-b2group['VMAG'], yerr = b2group['Vmagerr'], fmt = 'D', elinewidth  = 2, c = 'k', label = 'V')
plt.plot(x, func(x,popt[0],popt[1]), linestyle = '-', label = "[V]k' = %.3f"%popt[0], c = 'k',linewidth = 1)
plt.xlabel('Airmass', fontsize = 15)
k1V =  popt[0]
k1Verr = np.sqrt(np.diag(pcov))
popt, pcov = curve_fit(func, airlist, Bmaglist, sigma = ab2group['Bmagerr'])

plt.errorbar(x = a2group['airmass'],y = a2group['Bmag']-a2group['BMAG'], yerr = a2group['Bmagerr'], fmt = 'X', elinewidth  = 2, c = 'k')
plt.errorbar(x = b2group['airmass'],y = b2group['Bmag']-b2group['BMAG'], yerr = b2group['Bmagerr'], fmt = 'X', elinewidth  = 2, c = 'k', label  ='B')
plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[B]k' = %.3f"%popt[0], c = 'k', linewidth = 1)
k1B = popt[0]
k1Berr = np.sqrt(np.diag(pcov))
plt.ylabel(r'$mag_{obs} - mag_{ref}$', fontsize = 15)
plt.legend()
plt.show()



#%% k''
def main(targetlist):
    from astropy.table import vstack
    kVlist = []
    kBlist = []
    for target1 in targetlist:
        group_obj = tbl.group_by('object')
        a1group = tbl[tbl['object'] == target1]
        tbl1_ref = tbl_ref[tbl_ref['object'] == target1] 
        a1group['VMAG'] = tbl1_ref['VMAG']
        a1group['BMAG'] = tbl1_ref['B-V']+tbl1_ref['VMAG']
        a1group['B-V'] = tbl1_ref['B-V']
        
        
        air1list = a1group['airmass']
        V1maglist = a1group['Vmag']-a1group['VMAG']
        B1maglist = a1group['Bmag']-a1group['BMAG']
        
        from scipy.optimize import curve_fit
        
        def func(x, a, b):
            y = a * x + b
            return y
        
        def fit_show(func, tbl, title = None):
            airlist = tbl['airmass']
            Vmaglist = tbl['Vmag']-tbl['VMAG']
            Bmaglist = tbl['Bmag']-tbl['BMAG']
        
            x = np.linspace(np.min(airlist),np.max(airlist))
            popt, pcov = curve_fit(func, airlist, Vmaglist)
            plt.figure(figsize  =(6,4), dpi = 300)
            plt.title(f'STD-{title}')
            plt.xlabel('airmass', fontsize = 15)
            plt.errorbar(x = tbl['airmass'],y = Vmaglist, yerr = tbl['Vmagerr'], c = 'k', fmt = 'D', elinewidth  = 3, label ='V')
            plt.plot(x, func(x,popt[0],popt[1]), linestyle = '-', c = 'k', label = "[V]k = %.2f"%popt[0])
            kV =  popt[0]
            plt.ylabel(r'$mag_{obs}-mag_{ref}$', fontsize  =15)
            popt, pcov = curve_fit(func, airlist, Bmaglist)
            plt.errorbar(x = tbl['airmass'],y = Bmaglist, yerr = tbl['Bmagerr'], fmt = 'X', c= 'k', elinewidth  = 3, label = 'B')
            plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', c= 'k', label = "[B]k = %.2f"%popt[0])
            kB = popt[0]
            plt.legend()
            
            return kV, kB
            
        kV, kB = fit_show(func, a1group, title = target1)
        kVlist.append(kV)
        kBlist.append(kB)
        plt.show()
        
    return kVlist, kBlist

#%% k''
kV, kB = main(['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5'])
AkV, AkB = main(['a1','a2','a3','a4','a5'])
BkV, BkB = main(['b1','b2','b3','b4','b5'])

A_tbl_ref = tbl_ref[:5]
B_tbl_ref = tbl_ref[5:]
Ak1B = AkB[1]
Ak1V = AkV[1]
Bk1B = BkB[1]
Bk1V = BkV[1]

#k'' estimation from different color stars in the same field
Ak2B = (AkB[4]-AkB[0])/(A_tbl_ref[4]['B-V']-A_tbl_ref[0]['B-V'])
Ak2V = (AkV[4]-AkV[0])/(A_tbl_ref[4]['B-V']-A_tbl_ref[0]['B-V'])
Bk2B = (BkB[4]-BkB[0])/(B_tbl_ref[4]['B-V']-B_tbl_ref[0]['B-V'])
Bk2V = (BkV[4]-BkV[0])/(B_tbl_ref[4]['B-V']-B_tbl_ref[0]['B-V'])
#%%
def show(reftbl, result, marker = 'o'):
    plt.figure(figsize = (6,4), dpi = 500)
    x = np.linspace(np.min(reftbl['B-V']),np.max(reftbl['B-V']))
    popt, pcov = curve_fit(func, reftbl['B-V'], result)
    plt.scatter(reftbl['B-V'], result, c= 'k', marker = marker, label = 'Field A')
    #plt.scatter(reftbl['B-V'][:5], result[:5], c= 'k', marker = marker, label = 'Field A')
    #plt.scatter(reftbl['B-V'][5:], result[5:], c= 'r', marker = marker, label = 'Field B')
    plt.xlabel('B-V', fontsize = 15)
    plt.ylabel(r'$K_B$', fontsize = 15)
    plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[B]k'' = %.3f"%popt[0], c = 'k', linewidth  =1)
    err = np.sqrt(np.diag(pcov))
    plt.legend()
    plt.show()
    return popt[0]
#%%
Ak2B = show(A_tbl_ref, AkB, 'X')
Ak2V = show(A_tbl_ref, AkV)
Bk2B = show(B_tbl_ref, BkB, 'X')
Bk2V = show(B_tbl_ref, BkV)
k2B = show(tbl_ref, kB, 'X')
k2V = show(tbl_ref, kV, 'D')

#%% k''
plt.figure(figsize = (6,4), dpi = 500)
x = np.linspace(np.min(tbl_ref['B-V']),np.max(tbl_ref['B-V']))
popt, pcov = curve_fit(func, tbl_ref['B-V'], kV)
#plt.scatter(reftbl['B-V'], result, c= 'r', marker = marker, label = 'Field B')
plt.plot(x, func(x,popt[0],popt[1]), linestyle = '-', label = "[V]k'' = %.3f"%popt[0], c = 'k', linewidth  =1)


plt.xlabel('B-V', fontsize = 15)
plt.ylabel(r'$K$', fontsize = 15)
popt, pcov = curve_fit(func, tbl_ref['B-V'], kB)
#plt.scatter(reftbl['B-V'], result, c= 'r', marker = marker, label = 'Field B')


plt.xlabel('B-V', fontsize = 15)
plt.ylabel(r'$K$', fontsize = 15)
plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[B]k'' = %.3f"%popt[0], c = 'k', linewidth  =1)
plt.scatter(tbl_ref[5:]['B-V'], kV[5:], c= 'r', marker = 'D')
plt.scatter(tbl_ref[:5]['B-V'], kV[:5], c= 'k', marker = 'o', label = 'Field A')
plt.scatter(tbl_ref[5:]['B-V'], kV[5:], c= 'r', marker = 'o', label = 'Field B')
plt.scatter(tbl_ref[:5]['B-V'], kV[:5], c= 'k', marker = 'D', label = 'V')
plt.scatter(tbl_ref[:5]['B-V'], kB[:5], c= 'k', marker = 'X', label = 'B')
plt.scatter(tbl_ref[5:]['B-V'], kB[5:], c= 'r', marker = 'X')
plt.legend()
plt.show()
#%% a
from astropy.table import join
tbl['Vmag']

a_idx = [tbl['object'][i].startswith('a') for i in range(len(tbl))]
b_idx = [tbl['object'][i].startswith('b') for i in range(len(tbl))]

a_tbl = tbl[a_idx]
b_tbl = tbl[b_idx]

tbl['cor_Vmag'] = tbl['Vmag'] - k1V *tbl['airmass'] - k2V *tbl['BVmag']*tbl['airmass']
tbl['cor_Bmag'] = tbl['Bmag'] - k1B *tbl['airmass'] - k2B *tbl['BVmag']*tbl['airmass']
#a_tbl['cor_Vmag'] = a_tbl['Vmag'] - Ak1V *a_tbl['airmass'] - Ak2V *a_tbl['BVmag']*a_tbl['airmass']
##a_tbl['cor_Bmag'] = a_tbl['Bmag'] - Ak1V *a_tbl['airmass'] - Ak2V *a_tbl['BVmag']*a_tbl['airmass']
#b_tbl['cor_Vmag'] = b_tbl['Vmag'] - Bk1V *b_tbl['airmass'] - Bk2V *b_tbl['BVmag']*b_tbl['airmass']
#b_tbl['cor_Bmag'] = b_tbl['Bmag'] - Bk1V *b_tbl['airmass'] - Bk2V *b_tbl['BVmag']*b_tbl['airmass']
#a_tbl = a_tbl[a_tbl['airmass'] == 1.0]
#b_tbl = b_tbl[b_tbl['airmass'] == 1.5]

#%%
def show1(reftbl, result, marker = '.'):
    x = np.linspace(np.min(reftbl['B-V']),np.max(reftbl['B-V']))
    popt, pcov = curve_fit(func, reftbl['B-V'], result)
    plt.figure(figsize = (6,4), dpi = 300)
    plt.scatter(reftbl['B-V'], result, c= 'k', marker = marker)
    plt.xlabel('B-V', fontsize = 15)
    plt.ylabel(r'$magV_{cor}-magV_{ref}$', fontsize = 15)
    plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[V]a = %.2f"%popt[0], c = 'k')
    plt.legend()
    plt.show()
    return popt[0]

#AaB = show1(a_tbl, a_tbl['cor_Bmag'] - a_tbl['BMAG'])
#AaV = show1(a_tbl, a_tbl['cor_Vmag'] - a_tbl['VMAG'])
#BaB = show1(b_tbl, b_tbl['cor_Bmag'] - b_tbl['BMAG'])
#BaV = show1(b_tbl, b_tbl['cor_Vmag'] - b_tbl['VMAG'])
aB = show1(tbl, tbl['cor_Bmag']-tbl['BMAG'], 'X')
aV = show1(tbl, tbl['cor_Vmag'] - tbl['VMAG'], 'D')
#%% c
def show1(reftbl, result, marker = '.'):
    x = np.linspace(np.min(reftbl['B-V']),np.max(reftbl['B-V']))
    popt, pcov = curve_fit(func, reftbl['B-V'], result)
    plt.figure(figsize = (6,4), dpi = 300)
    plt.scatter(reftbl['B-V'], result, c= 'k', marker = marker)
    plt.xlabel('B-V', fontsize = 15)
    plt.ylabel(r'$(b-v)_{cor}$', fontsize = 15)
    plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "c = %.3f"%popt[0], c = 'k')
    plt.legend()
    plt.show()
    return popt[0]

#AaB = show1(a_tbl, a_tbl['cor_Bmag'] - a_tbl['BMAG'])
#AaV = show1(a_tbl, a_tbl['cor_Vmag'] - a_tbl['VMAG'])
#BaB = show1(b_tbl, b_tbl['cor_Bmag'] - b_tbl['BMAG'])
#BaV = show1(b_tbl, b_tbl['cor_Vmag'] - b_tbl['VMAG'])
aB = show1(tbl, tbl['cor_Bmag']-tbl['BMAG'], 'X')
aV = show1(tbl, tbl['cor_Vmag'] - tbl['VMAG'], 'D')
k1bv = k1B-k1V
k2bv = k2B-k2V
tbl['cor_BVmag'] = tbl['cor_Bmag']-tbl['cor_Vmag']
cor_BVmag = tbl['BVmag']-k1bv*tbl['airmass']-k2bv*tbl['airmass']*tbl['B-V']
cor_BVmag
tbl['cor_BVmag']
c, pcov = show1(tbl, tbl['cor_BVmag'])
np.sqrt(np.diag(pcov))
c
plt.scatter(tbl['B-V'], tbl['cor_BVmag'])
plt.scatter(tbl['B-V'], cor_BVmag)
np.median(tbl['cor_BVmag']-tbl['B-V'])
#%%
"""
all_tbl = tbl
plt.figure(figsize = (6,4), dpi = 500)

airmasslist = all_tbl.group_by('airmass')
for tbl in airmasslist.groups:
    #tbl = tbl[tbl['airmass'] == 2.0]
    x = np.linspace(np.min(tbl_ref['B-V']),np.max(tbl_ref['B-V']))
    popt, pcov = curve_fit(func, tbl['B-V'], tbl['cor_Vmag'] - tbl['VMAG'], sigma = tbl['Vmagerr'])
    #plt.scatter(reftbl['B-V'], result, c= 'r', marker = marker, label = 'Field B')
    aV = popt[0]
    aVerr = np.sqrt(np.diag(pcov))
    plt.scatter(tbl['B-V'], tbl['cor_Vmag'] - tbl['VMAG'], marker = 'D')
    plt.xlabel('B-V', fontsize = 15)
    plt.ylabel(r'$mag_{cor} - mag_{ref}$', fontsize = 15)
    plt.plot(x, func(x,popt[0],popt[1]), linestyle = '-', label = "[V]a = %.3f"%popt[0], linewidth  =1)
for tbl in airmasslist.groups:
    popt, pcov = curve_fit(func, tbl['B-V'], tbl['cor_Bmag']-tbl['BMAG'], sigma = tbl['Bmagerr'])
    #plt.scatter(reftbl['B-V'], result, c= 'r', marker = marker, label = 'Field B')
    aB = popt[0]
    aBerr = np.sqrt(np.diag(pcov))
    
    plt.scatter(tbl['B-V'], tbl['cor_Bmag']-tbl['BMAG'], marker = 'X')
    
    plt.xlabel('B-V', fontsize = 15)
    plt.ylabel(r'$mag_{cor} - mag_{ref}$', fontsize = 15)
    plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[B]a = %.3f"%popt[0], linewidth  =1)
all_tbl.group_by('airmass')
airmasslist.groups
plt.legend()
plt.show()
aBerr
aVerr
"""
#%% ZP








#%%
# Defining functions (for convenience)
# importing necessary modules
import numpy as np
import glob, os
import pandas as pd
from sklearn import linear_model
import statsmodels.api as sm
from matplotlib import pyplot as plt
# Plot - Observed values vs. Fitted values
def plot_comparison(input_data, fitted_data):
    arr0 = np.linspace(-5.0, 0.0, 1000)
    min_limit = np.minimum(input_data.min(), fitted_data.min()) - 0.2
    max_limit = np.maximum(input_data.max(), fitted_data.max()) + 0.2

    fig, ax = plt.subplots(figsize=(5,5))
    ax.plot(arr0, arr0, 'r--', linewidth=1.5, alpha=0.6)
    ax.plot(input_data, fitted_data, 'o', color='blue', ms=4.0)
    ax.tick_params(axis='both', labelsize=12.0)
    ax.set_xlabel(r"Observed $r-R$", fontsize=12.0)
    ax.set_ylabel(r"Fitted $r-R$", fontsize=12.0)
    ax.set_xlim([min_limit, max_limit])
    ax.set_ylim([min_limit, max_limit])
    plt.tight_layout()

# Plot - Observed values vs. Residuals
def plot_residuals(input_data, residuals):
    arr0 = np.linspace(-5.0, 0.0, 1000)
    min_limit = input_data.min() - 0.2
    max_limit = input_data.max() + 0.2
    RMSE = np.sqrt(np.sum(residuals**2) / len(input_data))

    fig, ax = plt.subplots(figsize=(5,5))
    ax.plot(arr0, np.zeros_like(arr0), 'r--', linewidth=1.5, alpha=0.6)
    ax.plot(input_data, residuals, 'o', color='blue', ms=4.0)
    ax.tick_params(axis='both', labelsize=12.0)
    ax.set_xlabel(r"Observed $r-R$", fontsize=12.0)
    ax.set_ylabel("Residuals", fontsize=12.0)
    ax.set_xlim([min_limit, max_limit])
    ax.set_ylim([-1.5, 1.5])
    ax.text(0.05, 0.95, f"RMS Error = {RMSE:.2f}", fontsize=13.0, fontweight='bold',
            transform=ax.transAxes, ha='left', va='top')
    plt.tight_layout()

# Printing the summary of model
def summary_model(x, y, e_y):
    Xm = sm.add_constant(x)
    model = sm.WLS(y.astype('float'), Xm.astype('float'), weights=1/e_y**2).fit() 
    print_model = model.summary()
    print(print_model)
    
#%%
from sklearn.linear_model import LinearRegression

import pandas as pd

tbl_fit = tbl
tbl_pd = tbl_fit.to_pandas()
X = tbl_pd[['airmass','B-V']]
YB = tbl_pd['Bmag'] - tbl_pd['BMAG']
e_YB = tbl_pd['Bmagerr']
YV = tbl_pd['Vmag'] - tbl_pd['VMAG']
e_YV = tbl_pd['Vmagerr']

#%%
#B
mlr = LinearRegression()
mlr.fit(X,YB, 1/e_YB**2)
fitted_YB = mlr.predict(X)
resiB = YB-fitted_YB
summary_model(X, YB, e_YB)
plot_comparison(fitted_YB, YB)
plot_residuals(YB, resiB)
fit_kB, fit_aB = mlr.coef_
fit_zpB = mlr.intercept_
#%%
#V
mlr = LinearRegression()
mlr.fit(X,YV, 1/e_YV**2)
fitted_YV = mlr.predict(X)
resiV = YV-fitted_YV
summary_model(X, YV, e_YV)
plot_comparison(fitted_YV, YV)
plot_residuals(YV, resiV)
fit_kV, fit_aV = mlr.coef_
fit_zpV = mlr.intercept_

#%%
tbl_fit = tbl
tbl_pd = tbl_fit.to_pandas()
tbl_pd['(B-V)X'] = tbl_pd['airmass']*tbl_pd['B-V']
X = tbl_pd[['airmass','B-V','(B-V)X']]
YB = tbl_pd['Bmag'] - tbl_pd['BMAG']
e_YB = tbl_pd['Bmagerr']
YV = tbl_pd['Vmag'] - tbl_pd['VMAG']
e_YV = tbl_pd['Vmagerr']

mlr = LinearRegression()
mlr.fit(X,YB, 1/e_YB**2)
fitted_YB = mlr.predict(X)
resiB = YB-fitted_YB
summary_model(X, YB, e_YB)
plot_comparison(fitted_YB, YB)
plot_residuals(YB, resiB)
fit_kB, fit_aB = mlr.coef_
fit_zpB = mlr.inter0.cept_
#%%
mlr = LinearRegression()
mlr.fit(X,YV, 1/e_YV**2)
fitted_YV = mlr.predict(X)
resiV = YV-fitted_YV
summary_model(X, YV, e_YV)
plot_comparison(fitted_YV, YV)
plot_residuals(YV, resiV)
fit_kV, fit_aV = mlr.coef_
fit_zpV = mlr.intercept_

#%%
tbl_fit = tbl
tbl_pd = tbl_fit.to_pandas()
tbl_pd['(B-V)X'] = tbl_pd['airmass']*tbl_pd['B-V']
X = tbl_pd[['airmass','B-V','(B-V)X']]
YB = tbl_pd['BVmag']
e_YB = tbl_pd['Bmagerr']

mlr = LinearRegression()
mlr.fit(X,YB, 1/e_YB**2)
fitted_YB = mlr.predict(X)
resiB = YB-fitted_YB
summary_model(X, YB, e_YB)
plot_comparison(fitted_YB, YB)
plot_residuals(YB, resiB)
fit_kB, fit_aB = mlr.coef_
fit_zpB = mlr.inter0.cept_


#%%
tbl['Vmag']
def result(tbl, filter_, a, k1, k2):
    mag = tbl[f'{filter_}mag'] -( a*tbl['B-V'] +k1 *tbl['airmass'] + k2 *tbl['BVmag']*tbl['airmass'])
    return mag
def result2(tbl, a, k, zp):
    mag = tbl['Vmag'] -( a*tbl['B-V'] +k *tbl['airmass'] ) - zp
    return mag
def result3(tbl, filter_, a, k1, k2, zp):
    mag = tbl[f'{filter_}mag'] -( a*tbl['B-V'] +k1 *tbl['airmass'] + k2 *tbl['BVmag']*tbl['airmass']) + zp
    return mag
tbl['real_Vmag'] = result(tbl, 'V', aV, k1V, k2V)
tbl['real_Bmag'] = result(tbl, 'B', aB, k1B, k2B)
#%%
def func3(x, a):
    y = x + a
    return y
popt, pcov= curve_fit(func3, tbl['VMAG'], tbl['real_Vmag'], sigma = tbl['Vmagerr'])
zpVerr = np.sqrt(np.diag(pcov))
zpVerr = np.std(tbl['real_Vmag']-tbl['VMAG'])
zpV = np.median(tbl['real_Vmag']-tbl['VMAG'])
plt.figure(figsize = (6,4), dpi = 500)
plt.errorbar(tbl['Vmag']-tbl['VMAG'], tbl['real_Vmag']-tbl['VMAG'], yerr = tbl['Vmagerr'], fmt= 'D', capsize  =3, c = 'k', label  ='V')
plt.ylim(-0.1,-0.4)
plt.xlabel(r'$mag_{obs} - mag_{ref}$', fontsize = 15)
plt.ylabel(r'$mag_{cor} - mag_{ref}$', fontsize = 15)
plt.legend()
plt.text(-0.18,-0.37,'ZP = %.3f'%zpV, fontsize = 12)
zpV
zpVerr
#%%
popt, pcov= curve_fit(func3, tbl['BMAG'], tbl['real_Bmag'], sigma  =tbl['Bmagerr'])
zpB = np.median(tbl['real_Bmag']-tbl['BMAG'])
zpBerr = np.sqrt(np.diag(pcov))
zpBerr = np.std(tbl['real_Bmag']-tbl['BMAG'])
plt.figure(figsize = (6,4), dpi = 500)
plt.errorbar(tbl['BMAG']-tbl['Bmag'], tbl['real_Bmag']-tbl['BMAG'], yerr=  tbl['Bmagerr'], fmt= 'X', c = 'k', capsize =3, label = 'B')
plt.ylim(0.3, -0.2)
plt.xlabel(r'$mag_{obs} - mag_{ref}$', fontsize = 15)
plt.ylabel(r'$mag_{cor} - mag_{ref}$', fontsize = 15)
plt.text(-0.65,-0.15,'ZP = %.3f'%zpB, fontsize = 12)
plt.legend()
tbl['Bmagerr']
zpB
zpBerr
#%%
tbl[(tbl['BMAG']-tbl['real_Bmag'] < -0.1)]['airmass']
k2V
tbl['real_Vmag'] = result3(tbl, 'V', aV, k1V, k2V, zpV)
tbl['real_Bmag'] = result3(tbl, 'B', aB, k1B, k2B, zpB)
fit_result = result2(tbl, fit_a, fit_k, fit_zp)
#plt.scatter(tbl['B-V'], tbl['VMAG']-tbl['real_Vmag'])
#plt.scatter(tbl['B-V'], tbl['VMAG']-tbl['Vmag'])
plt.scatter(tbl['VMAG']-tbl['Vmag'], tbl['VMAG']-tbl['real_Vmag'])
plt.scatter(tbl['BMAG']-tbl['Bmag'], tbl['BMAG']-tbl['real_Bmag'])
plt.scatter(tbl['VMAG']-tbl['Vmag'],tbl['VMAG']-fit_result)
plt.xlim(-0.1,0.5)
plt.ylim(-0.1,0.5)


tbl
#%%
def result4(tbl, k1bv, k2bv):
    mag = tbl[f'B-V'] - k1bv*tbl['airmass'] -k2bv*tbl['B-V']*tbl['airmass']
    return mag
k1bv = k1B-k1V
k2bv = k2B-k2V
tbl['cor_BVmag'] = result4(tbl, k1bv, k2bv)
tbl_group = tbl.group_by('airmass')
plt.scatter(tbl['B-V'],tbl['cor_BVmag']-tbl['B-V'])
for table in tbl_group.groups:
    plt.errorbar(table['BVmag'],table['B-V']-table['cor_BVmag'], yerr = np.sqrt(table['Bmagerr']**2+table['Vmagerr']**2), fmt= 'o')
    x = np.linspace(np.min(table['BVmag']),np.max(table['BVmag']))
    popt, pcov = curve_fit(func, table['BVmag'], table['B-V']-table['cor_BVmag'], sigma = np.sqrt(table['Bmagerr']**2+table['Vmagerr']**2))
    #plt.scatter(reftbl['B-V'], result, c= 'r', marker = marker, label = 'Field B')
    aBV = popt[0]
    aBVerr = np.sqrt(np.diag(pcov))

    plt.xlabel('B-V', fontsize = 15)
    plt.ylabel(r'$mag_{cor} - mag_{ref}$', fontsize = 15)
    plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[B]a = %.3f"%popt[0], linewidth  =1)






#%% k 
plt.figure(figsize = (6,4), dpi = 500)
x = np.linspace(np.min(tbl['airmass']),np.max(tbl['airmass']))
popt, pcov = curve_fit(func, tbl['airmass'], tbl['Vmag']-tbl['VMAG'], sigma = tbl['Vmagerr'])
plt.errorbar(x = tbl['airmass'],y = tbl['Vmag']-tbl['VMAG'], yerr = tbl['Vmagerr'], fmt = 'D', elinewidth  = 2, c = 'k')
plt.plot(x, func(x,popt[0],popt[1]), linestyle = '-', label = "[V]k = %.3f"%popt[0], c = 'k',linewidth = 1)
plt.xlabel('Airmass', fontsize = 15)
k1V =  popt[0]
k1Verr = np.sqrt(np.diag(pcov))

popt, pcov = curve_fit(func, tbl['airmass'], tbl['Bmag']-tbl['BMAG'], sigma = tbl['Bmagerr'])
plt.errorbar(x = tbl['airmass'],y = tbl['Bmag']-tbl['BMAG'], yerr = tbl['Bmagerr'], fmt = 'X', elinewidth  = 2, c = 'k')
plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[B]k = %.3f"%popt[0], c = 'k', linewidth = 1)
k1B = popt[0]
k1Berr = np.sqrt(np.diag(pcov))
plt.ylabel(r'$mag_{obs} - mag_{ref}$', fontsize = 15)
plt.legend()
plt.show()

#%% a
from astropy.table import join
a_idx = [tbl['object'][i].startswith('a') for i in range(len(tbl))]
b_idx = [tbl['object'][i].startswith('b') for i in range(len(tbl))]
a_tbl = tbl[a_idx]
b_tbl = tbl[b_idx]

tbl['cor_Vmag'] = tbl['Vmag'] - k1V *tbl['airmass'] 
tbl['cor_Bmag'] = tbl['Bmag'] - k1B *tbl['airmass'] 

plt.figure(figsize = (6,4), dpi = 500)
x = np.linspace(np.min(tbl_ref['B-V']),np.max(tbl_ref['B-V']))
popt, pcov = curve_fit(func, tbl['B-V'], tbl['cor_Vmag'] - tbl['VMAG'], sigma = tbl['Vmagerr'])
aV = popt[0]
aVerr = np.sqrt(np.diag(pcov))
plt.scatter(tbl['B-V'], tbl['cor_Vmag'] - tbl['VMAG'], c= 'k', marker = 'D', label = 'V')
plt.xlabel('B-V', fontsize = 15)
plt.ylabel(r'$mag_{cor} - mag_{ref}$', fontsize = 15)
plt.plot(x, func(x,popt[0],popt[1]), linestyle = '-', label = "[V]a = %.3f"%popt[0], c = 'k', linewidth  =1)
#%%
popt, pcov = curve_fit(func, tbl['B-V'], tbl['cor_Bmag']-tbl['BMAG'], sigma = tbl['Bmagerr'])
aB = popt[0]
aBerr = np.sqrt(np.diag(pcov))
plt.scatter(tbl['B-V'], tbl['cor_Bmag']-tbl['BMAG'], c= 'k', marker = 'X', label = 'B')
plt.xlabel('B-V', fontsize = 15)
plt.ylabel(r'$mag_{cor} - mag_{ref}$', fontsize = 15)
plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[B]a = %.3f"%popt[0], c = 'k', linewidth  =1)

plt.legend()
plt.show()


#%% ZP
def result(tbl, filter_, a, k1, k2):
    mag = tbl[f'{filter_}mag'] -( a*tbl['B-V'] +k1 *tbl['airmass'] + k2 *tbl['BVmag']*tbl['airmass'])
    return mag
def result2(tbl, filter_, a, k):
    mag = tbl[f'{filter_}mag'] -( a*tbl['B-V'] +k *tbl['airmass'] )
    return mag
def result3(tbl, filter_, a, k1, k2, zp):
    mag = tbl[f'{filter_}mag'] -( a*tbl['B-V'] +k1 *tbl['airmass'] + k2 *tbl['BVmag']*tbl['airmass']) + zp
    return mag
tbl['real_Vmag'] = result2(tbl, 'V', aV, k1V)
tbl['real_Bmag'] = result2(tbl, 'B', aB, k1B)

def func3(x, a):
    y = x + a
    return y
popt, pcov= curve_fit(func3, tbl['real_Vmag'], tbl['VMAG'], sigma = tbl['Vmagerr'])
#zpVerr = np.sqrt(np.diag(pcov))
zpVerr = np.std(tbl['real_Vmag']-tbl['VMAG'])
zpV = np.median(tbl['real_Vmag']-tbl['VMAG'])
plt.figure(figsize = (6,4), dpi = 500)
plt.errorbar(tbl['Vmag']-tbl['VMAG'], tbl['real_Vmag']-tbl['VMAG'], yerr = tbl['Vmagerr'], fmt= 'D', capsize  =3, c = 'k', label  ='V')
plt.ylim(-0.1,-0.4)
plt.xlabel(r'$mag_{obs} - mag_{ref}$', fontsize = 15)
plt.ylabel(r'$mag_{cor} - mag_{ref}$', fontsize = 15)
plt.legend()
plt.text(-0.18,-0.37,'ZP = %.3f'%zpV, fontsize = 12)
zpV
#%%
popt, pcov= curve_fit(func3, tbl['BMAG'], tbl['real_Bmag'], sigma  =tbl['Bmagerr'])
tbl['real_Bmag']
tbl['BMAG']
zpB = np.median(tbl['real_Bmag']-tbl['BMAG'])
zpBerr = np.std(tbl['real_Bmag']-tbl['BMAG'])
#zpBerr = np.sqrt(np.diag(pcov))
plt.figure(figsize = (6,4), dpi = 500)
plt.errorbar(tbl['BMAG']-tbl['Bmag'], tbl['real_Bmag']-tbl['BMAG'], yerr=  tbl['Bmagerr'], fmt= 'X', c = 'k', capsize =3, label = 'B')
plt.ylim(0.3, -0.2)
plt.xlabel(r'$mag_{obs} - mag_{ref}$', fontsize = 15)
plt.ylabel(r'$mag_{cor} - mag_{ref}$', fontsize = 15)
plt.text(-0.65,-0.15,'ZP = %.3f'%zpB, fontsize = 12)
plt.legend()
zpBerr

#%% c
def show1(reftbl, result, marker = '.'):
    x = np.linspace(np.min(reftbl['B-V']),np.max(reftbl['B-V']))
    popt, pcov = curve_fit(func, reftbl['B-V'], result, sigma = np.sqrt(reftbl['Bmagerr']**2+reftbl['Vmagerr']**2))
    plt.figure(figsize = (6,4), dpi = 300)
    plt.scatter(reftbl['B-V'], result, c= 'k', marker = marker)
    plt.xlabel('B-V', fontsize = 15)
    plt.ylabel(r'$(b-v)_{cor}$', fontsize = 15)
    plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "c = %.3f"%popt[0], c = 'k')
    plt.legend()
    plt.show()
    return popt[0], pcov
#%%

tbl['cor_Vmag'] = tbl['Vmag'] - k1V *tbl['airmass'] 
tbl['cor_Bmag'] = tbl['Bmag'] - k1B *tbl['airmass'] 
#AaB = show1(a_tbl, a_tbl['cor_Bmag'] - a_tbl['BMAG'])
#AaV = show1(a_tbl, a_tbl['cor_Vmag'] - a_tbl['VMAG'])
#BaB = show1(b_tbl, b_tbl['cor_Bmag'] - b_tbl['BMAG'])
#BaV = show1(b_tbl, b_tbl['cor_Vmag'] - b_tbl['VMAG'])
aB = show1(tbl, tbl['cor_Bmag']-tbl['BMAG'], 'X')
aV = show1(tbl, tbl['cor_Vmag'] - tbl['VMAG'], 'D')
k1bv = k1B-k1V
k2bv = k2B-k2V
tbl['cor_BVmag'] = tbl['cor_Bmag']-tbl['cor_Vmag']
cor_BVmag = tbl['BVmag']-k1bv*tbl['airmass']
cor_BVmag
c, pcov = show1(tbl, tbl['cor_BVmag'])
np.sqrt(np.diag(pcov))
plt.scatter(tbl['B-V'], tbl['cor_BVmag'])
plt.scatter(tbl['B-V'], cor_BVmag)
plt.scatter(tbl['B-V'],tbl['cor_BVmag']-c*tbl['B-V'])
np.median(tbl['cor_BVmag']-c*tbl['B-V'])
np.std(tbl['cor_BVmag']-c*tbl['B-V'])



#%%
ik1v, ie_k1v = 0.130, 0.005
ik2v, ie_k2v = 0.019, 0.004
iav, ie_av = -0.054, 0.009
izpv, ie_zpv = -0.248, 0.008
ik1b, ie_k1b = 0.242, 0.018
ik2b, ie_k2b = 0.028, 0.018
iab, ie_ab = -0.114, 0.033
izpb, ie_zpb = 0.062, 0.031
ic = 0.939
ik1bv = ik1b-ik1v
ie_k1bv = np.sqrt(ie_k1v**2+ie_k1b**2)
ie_k2bv = np.sqrt(ie_k2v**2+ie_k2b**2)

from scipy.optimize import fsolve
target = tbl[-4:]
iVmaglist = []
iBmaglist = []
iVmagerrlist = []
iBmagerrlist = []
#%%
for i in range(4):
    target['airmass']
    airmass = target[i]['airmass']
    bv = target[i]['BVmag']
    vobs = target[i]['Vmag']
    bobs = target[i]['Bmag']
    e_vobs = target[i]['Vmagerr']
    e_bobs = target[i]['Bmagerr']
    def equations(var):
        B, V = var
        f1 = izpb + ik1b*airmass + ik2b*airmass*bv + iab*bv - (bobs-B)
        f2 = izpv + ik1v*airmass + ik2v*airmass*bv + iav*bv - (vobs-V)
        return [f1, f2]
    
    solution, infodict, ier, mesg = fsolve(equations, (bobs,vobs), full_output = True)
    
    BVcolor = solution[0]-solution[1]
    dBdzero = -1./(iab+1)
    dBdk = -airmass/(iab+1)
    dBdc = -BVcolor
    dBdv = 1./(iab+1)
    
    B_err = np.sqrt((dBdzero*ie_zpb)**2 +
                    (dBdk*ie_k1b)**2 +
                    (dBdc*ie_ab)**2 +
                    (dBdv*e_bobs)**2)
    
    BVcolor = solution[0]-solution[1]
    dVdzero = -1./(iav+1)
    dVdk = -airmass/(iav+1)
    dVdc = -BVcolor
    dVdv = 1./(iav+1)
    
    V_err = np.sqrt((dVdzero*ie_zpv)**2 +
                    (dVdk*ie_k1v)**2 +
                    (dVdc*ie_av)**2 +
                    (dVdv*e_vobs)**2)
    
    print(f'object = {target[i]["object"]}')
    print(f'Bmag = {target[i]["Bmag"]}')
    print(f'Bmagerr = {target[i]["Bmagerr"]}')
    print(f'Vmag = {target[i]["Vmag"]}')
    print(f'Vmagerr = {target[i]["Vmagerr"]}')
    print(f'BMAG = {solution[0]}')
    print(f'BMAGerr = {B_err}')
    print(f'VMAG = {solution[1]}')
    print(f'VMAGerr = {V_err}')
    print(f'Standard B-V = {BVcolor}')
    print(f'Standard B-Verr = {np.sqrt(B_err**2+V_err**2)}')
    iBmaglist.append(solution[0])
    iBmagerrlist.append(B_err)
    iVmaglist.append(solution[1])
    iVmagerrlist.append(V_err)
#%% Check!




#%%
k1v, e_k1v = 0.128, 0.001
k2v, e_k2v = 0.018, 0.001
av, e_av = -0.047, 0.003
zpv, e_zpv = -0.260, 0.010
k1b, e_k1b = 0.247, 0.034
k2b, e_k2b = 0.031, 0.023
ab, e_ab = -0.125, 0.026
zpb, e_zpb = 0.033, 0.034
c = 0.933
k1bv = k1b-k1v
e_k1bv = np.sqrt(e_k1v**2+e_k1b**2)
e_k2bv = np.sqrt(e_k2v**2+e_k2b**2)

from scipy.optimize import fsolve
target = tbl[-4:]
Vmaglist = []
Bmaglist = []
Vmagerrlist = []
Bmagerrlist = []
#%%
for i in range(4):
    target['airmass']
    airmass = target[i]['airmass']
    bv = target[i]['BVmag']
    vobs = target[i]['Vmag']
    bobs = target[i]['Bmag']
    e_vobs = target[i]['Vmagerr']
    e_bobs = target[i]['Bmagerr']
    def equations(var):
        B, V = var
        f1 = zpb + k1b*airmass + k2b*airmass*bv + ab*bv - (bobs-B)
        f2 = zpv + k1v*airmass + k2v*airmass*bv + av*bv - (vobs-V)
        return [f1, f2]
    
    solution, infodict, ier, mesg = fsolve(equations, (bobs,vobs), full_output = True)
    solution
    vobs
    tbl_ref
    tbl[2]
    
    BVcolor = solution[0]-solution[1]
    dBdzero = -1./(ab+1)
    dBdk = -airmass/(ab+1)
    dBdc = -BVcolor
    dBdv = 1./(ab+1)
    
    B_err = np.sqrt((dBdzero*e_zpb)**2 +
                    (dBdk*e_k1b)**2 +
                    (dBdc*e_ab)**2 +
                    (dBdv*e_bobs)**2)
    
    BVcolor = solution[0]-solution[1]
    dVdzero = -1./(av+1)
    dVdk = -airmass/(av+1)
    dVdc = -BVcolor
    dVdv = 1./(av+1)
    
    V_err = np.sqrt((dVdzero*e_zpv)**2 +
                    (dVdk*e_k1v)**2 +
                    (dVdc*e_av)**2 +
                    (dVdv*e_vobs)**2)
    
    print(f'object = {target[i]["object"]}')
    print(f'Bmag = {target[i]["Bmag"]}')
    print(f'Bmagerr = {target[i]["Bmagerr"]}')
    print(f'Vmag = {target[i]["Vmag"]}')
    print(f'Vmagerr = {target[i]["Vmagerr"]}')
    print(f'BMAG = {solution[0]}')
    print(f'BMAGerr = {B_err}')
    print(f'VMAG = {solution[1]}')
    print(f'VMAGerr = {V_err}')
    print(f'Standard B-V = {BVcolor}')
    print(f'Standard B-Verr = {np.sqrt(B_err**2+V_err**2)}')
    
    Bmaglist.append(solution[0])
    Bmagerrlist.append(B_err)
    Vmaglist.append(solution[1])
    Vmagerrlist.append(V_err)
    Vmaglist
    
    
#%%

def check(mag_cor, bv_cor, airmass, a, k1, k2, zp):
    mag_obs = mag_cor + a*bv_cor + k1*bv_cor+k2*bv_cor*airmass + zp
    return mag_obs
check(solution[0], BVcolor, airmass, ab, k1b, k2b, zpb)
BVcolor
solution[0]
plt.scatter(Bmaglist, iBmaglist)


#%%
plt.figure(figsize = (6,4),  dpi = 300)
plt.errorbar(Bmaglist, iBmaglist, Bmagerrlist, iBmagerrlist, fmt = ',', capsize = 5, c = 'k')
popt, pcov = curve_fit(func, Bmaglist, iBmaglist, sigma = Bmagerrlist)
popt
xrange = np.linspace(np.min(Bmaglist), np.max(Bmaglist),100)
reduced_chi_2 = np.sum(((np.array(Bmaglist)-func(np.array(Bmaglist), popt[0], popt[1]))/Bmagerrlist)**2)/2
popt[1]
plt.plot(xrange, xrange, linestyle = '--', c='k', linewidth = 1, label = r'$\chi_{reduced}=%.2f$'%reduced_chi_2)
plt.xlabel(r'$mag_{cor, Graphical}$', fontsize = 15)
plt.ylabel(r'$mag_{cor, IRAF}$', fontsize = 15)
plt.legend()