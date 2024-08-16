
#%%
from astropy.io import ascii
import numpy as np
#%%
from astropy.table import Table
from gecko_rtsmaker import rtsmaker
#%%
tbl = ascii.read('./pointing_50Mpc_cbnuo.csv', format = 'csv')
tbl.sort('rate_sum')
tbl.reverse()
#%%
#%%
#start_i = 30
#end_i = 5419
#step_i = 10
i = 300
totlist = []
obslist = []
#for i in np.arange(30, 5400, 50):
cut_tbl = tbl[:i]
obj = cut_tbl['col0']
ra = cut_tbl['cen_RA'].round(5)
dec = cut_tbl['cen_Dec'].round(5)
rate_sum = cut_tbl['rate_sum']
glade_ID = cut_tbl['GLADE']

result_tbl = Table() 
result_tbl['obj'] = obj
result_tbl['ra'] = ra
result_tbl['dec'] = dec
result_tbl['priority'] = rate_sum
result_tbl['note'] = glade_ID
result_tbl['dist'] = 100
result_tbl['ebv'] = 0
result_tbl['muvmag'] = -20
result_tbl['maxaxis'] = 10
result_tbl['minaxis'] = 10
result_tbl.write('./gecko_target.dat', format = 'ascii')
rtsmaker(save_path = './', catpath  ='./gecko_target.dat', start = '2023/01/20', end = '2023/01/21', headname = 'GECKO', altlimit = 40, numlimit = 10000)
result = ascii.read('./GECKO-20230120-rts_vis-CBNUO.txt')
totlist.append(i)
obslist.append(len(result))

# %%
import matplotlib.pyplot as plt
# %%
plt.plot(totlist, obslist)
# %%
