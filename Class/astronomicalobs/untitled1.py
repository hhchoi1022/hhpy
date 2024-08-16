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
tbl_ref.write('tbl_ref.dat', format = 'ascii')
tbl = join(tbl_ref, tbl, keys = 'object')
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
    y = a * x + b
    return y
ab2group['Vmagerr']
x = np.linspace(np.min(airlist),np.max(airlist))
popt, pcov = curve_fit(func, airlist, Vmaglist)
plt.errorbar(x = ab2group['airmass'],y = Vmaglist, yerr = ab2group['Vmagerr'], fmt = '.', elinewidth  = 3)
plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[V]k' = %.2f"%popt[0])
k1V =  popt[0]
popt, pcov = curve_fit(func, airlist, Bmaglist)
plt.errorbar(x = ab2group['airmass'],y = Bmaglist, yerr = ab2group['Bmagerr'], fmt = '.', elinewidth  = 3)
plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[B]k' = %.2f"%popt[0])
k1B = popt[0]
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
            plt.title(title)
            plt.errorbar(x = tbl['airmass'],y = Vmaglist, yerr = tbl['Vmagerr'], fmt = '.', elinewidth  = 3)
            plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[V]k = %.2f"%popt[0])
            kV =  popt[0]
            popt, pcov = curve_fit(func, airlist, Bmaglist)
            plt.errorbar(x = tbl['airmass'],y = Bmaglist, yerr = tbl['Bmagerr'], fmt = '.', elinewidth  = 3)
            plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[B]k = %.2f"%popt[0])
            kB = popt[0]
            plt.legend()
            
            return kV, kB
            
        kV, kB = fit_show(func, a1group, title = target1)
        kVlist.append(kV)
        kBlist.append(kB)
        plt.show()
        
    return kVlist, kBlist

#%% k''
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

def show(reftbl, result):
    x = np.linspace(np.min(reftbl['B-V']),np.max(reftbl['B-V']))
    popt, pcov = curve_fit(func, reftbl['B-V'], result)
    plt.scatter(reftbl['B-V'], result)
    plt.xlabel('B-V')
    plt.ylabel('k')
    plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[V]k'' = %.2f"%popt[0], c = 'k')
    plt.legend()
    plt.show()
    return popt[0]
A_tbl_ref
Ak2B = show(A_tbl_ref, AkB)
Ak2V = show(A_tbl_ref, AkV)
Bk2B = show(B_tbl_ref, BkB)
Bk2V = show(B_tbl_ref, BkV)
k2B = np.mean([Ak2B,Bk2B])
k2V = np.mean([Ak2V,Bk2V])
Ak2B
AkB
BkB
k1B
k2B
#%% a
from astropy.table import join
tbl['Vmag']

a_idx = [tbl['object'][i].startswith('a') for i in range(len(tbl))]
b_idx = [tbl['object'][i].startswith('b') for i in range(len(tbl))]

a_tbl = tbl[a_idx]
b_tbl = tbl[b_idx]

tbl['cor_Vmag'] = tbl['Vmag'] - k1V *tbl['airmass'] - k2V *tbl['BVmag']*tbl['airmass']
tbl['cor_Bmag'] = tbl['Bmag'] - k1B *tbl['airmass'] - k2B *tbl['BVmag']*tbl['airmass']
a_tbl['cor_Vmag'] = a_tbl['Vmag'] - Ak1V *a_tbl['airmass'] - Ak2V *a_tbl['BVmag']*a_tbl['airmass']
a_tbl['cor_Bmag'] = a_tbl['Bmag'] - Ak1V *a_tbl['airmass'] - Ak2V *a_tbl['BVmag']*a_tbl['airmass']
b_tbl['cor_Vmag'] = b_tbl['Vmag'] - Bk1V *b_tbl['airmass'] - Bk2V *b_tbl['BVmag']*b_tbl['airmass']
b_tbl['cor_Bmag'] = b_tbl['Bmag'] - Bk1V *b_tbl['airmass'] - Bk2V *b_tbl['BVmag']*b_tbl['airmass']
a_tbl = a_tbl[a_tbl['airmass'] == 1.0]
b_tbl = b_tbl[b_tbl['airmass'] == 1.5]

plt.scatter(a_tbl['B-V'], a_tbl['VMAG']-a_tbl['cor_Vmag'])
plt.scatter(a_tbl['B-V'], a_tbl['BMAG']-a_tbl['cor_Bmag'])
plt.scatter(b_tbl['B-V'], b_tbl['VMAG']-b_tbl['cor_Vmag'])
plt.scatter(b_tbl['B-V'], b_tbl['BMAG']-b_tbl['cor_Bmag'])


#%%
def show1(reftbl, result):
    x = np.linspace(np.min(reftbl['B-V']),np.max(reftbl['B-V']))
    popt, pcov = curve_fit(func, reftbl['B-V'], result)
    plt.figure()
    plt.scatter(reftbl['B-V'], result)
    plt.xlabel('B-V')
    plt.ylabel('k')
    plt.plot(x, func(x,popt[0],popt[1]), linestyle = '--', label = "[V]a = %.2f"%popt[0], c = 'k')
    plt.legend()
    plt.show()
    return popt[0]

AaB = show1(a_tbl, a_tbl['cor_Bmag'] - a_tbl['BMAG'])
AaV = show1(a_tbl, a_tbl['cor_Vmag'] - a_tbl['VMAG'])
BaB = show1(b_tbl, b_tbl['cor_Bmag'] - b_tbl['BMAG'])
BaV = show1(b_tbl, b_tbl['cor_Vmag'] - b_tbl['VMAG'])
aB = show1(tbl, tbl['cor_Bmag']-tbl['BMAG'])
aV = show1(tbl, tbl['cor_Vmag'] - tbl['VMAG'])
aV
aB

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
mlr = LinearRegression()
tbl_fit
tbl_pd = tbl_fit.to_pandas()
X = tbl_pd[['airmass','B-V']]
Y = tbl_pd['Vmag'] - tbl_pd['VMAG']
tbl_pd
e_Y = tbl_pd['Bmagerr']
mlr.fit(X,Y, 1/e_Y**2)
fitted_Y = mlr.predict(X)
resi = Y-fitted_Y

summary_model(X, Y, e_Y)
plot_comparison(fitted_Y, Y)
plot_residuals(Y, resi)

fit_k, fit_a = mlr.coef_
fit_zp = mlr.intercept_
#%%
tbl['Vmag']
def result(tbl, a, k1, k2):
    mag = tbl['Vmag'] -( a*tbl['B-V'] +k1 *tbl['airmass'] + k2 *tbl['BVmag']*tbl['airmass'])
    return mag
def result2(tbl, a, k, zp):
    mag = tbl['Vmag'] -( a*tbl['B-V'] +k *tbl['airmass'] ) - zp
    return mag
def result3(tbl, a, k1, k2, zp):
    mag = tbl['Vmag'] -( a*tbl['B-V'] +k1 *tbl['airmass'] + k2 *tbl['BVmag']*tbl['airmass']) + zp
    return mag
tbl['real_Vmag'] = result(tbl, aV, k1V, k2V)
zp = np.median(tbl['VMAG']-tbl['real_Vmag'])
tbl['real_Vmag'] = result3(tbl, aV, k1V, k2V, zp)
fit_result = result2(tbl, fit_a, fit_k, fit_zp)
#plt.scatter(tbl['B-V'], tbl['VMAG']-tbl['real_Vmag'])
#plt.scatter(tbl['B-V'], tbl['VMAG']-tbl['Vmag'])
plt.scatter(tbl['VMAG']-tbl['Vmag'], tbl['VMAG']-tbl['real_Vmag'])
plt.scatter(tbl['VMAG']-tbl['Vmag'],tbl['VMAG']-fit_result)
plt.xlim(-0.1,0.5)
plt.ylim(-0.1,0.5)

