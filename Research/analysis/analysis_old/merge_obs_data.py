

#%%
import glob
from astropy.table import vstack
from astropy.io import ascii
from observedphot import ObservedPhot
from astropy.table import Table
#%%
HM236_key = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/*/*hostmwextin2.36.dat'
H236_key = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/*/*hostextin2.36.dat'
HM310_key = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/*/*hostmwextin3.10.dat'
H310_key = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/*/*hostextin3.10.dat'
Noextin_key = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/*/*noextin.dat'
MWextin_key = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/*/*mwextin3.10.dat'

HM236_filelist = glob.glob(HM236_key)
H236_filelist = glob.glob(H236_key)
HM310_filelist = glob.glob(HM310_key)
H310_filelist = glob.glob(H310_key)
Noextin_filelist = glob.glob(Noextin_key)
mwextin_filelist = glob.glob(MWextin_key)

savepath = '/Users/hhchoi1022/Gitrepo/Data/SN2021aefx/observation/lightcurve/'
#%%

HM236_tbl = Table()
for HM236 in HM236_filelist:
    tbl = ascii.read(HM236, format = 'fixed_width')
    HM236_tbl = vstack([HM236_tbl, tbl])
    HM236_tbl.write(savepath+'Alldata_H_M_cor236.dat', format = 'ascii.fixed_width', overwrite = True)

H236_tbl = Table()
for H236 in H236_filelist:
    tbl = ascii.read(H236, format = 'fixed_width')
    H236_tbl = vstack([H236_tbl, tbl])
    H236_tbl.write(savepath+'Alldata_H_cor236.dat', format = 'ascii.fixed_width', overwrite = True)

HM310_tbl = Table()
for HM310 in HM310_filelist:
    tbl = ascii.read(HM310, format = 'fixed_width')
    HM310_tbl = vstack([HM310_tbl, tbl])
    HM310_tbl.write(savepath+'Alldata_H_M_cor310.dat', format = 'ascii.fixed_width', overwrite = True)

H310_tbl = Table()
for H310 in H310_filelist:
    tbl = ascii.read(H310, format = 'fixed_width')
    H310_tbl = vstack([H310_tbl, tbl])
    H310_tbl.write(savepath+'Alldata_H_cor310.dat', format = 'ascii.fixed_width', overwrite = True)

Noextin_tbl = Table()
for Noextin in Noextin_filelist:
    tbl = ascii.read(Noextin, format = 'fixed_width')
    Noextin_tbl = vstack([Noextin_tbl, tbl])
    Noextin_tbl.write(savepath+'Alldata_No_cor.dat', format = 'ascii.fixed_width', overwrite = True)

MWextin_tbl = Table()
for Noextin in Noextin_filelist:
    tbl = ascii.read(Noextin, format = 'fixed_width')
    Noextin_tbl = vstack([Noextin_tbl, tbl])
    Noextin_tbl.write(savepath+'Alldata_M_cor310.dat', format = 'ascii.fixed_width', overwrite = True)

# %%
import matplotlib.pyplot as plt
plt.figure(figsize = (10,8))
data = ObservedPhot(data_tbl = HM310_tbl)
data.show_lightcurve(label = True)
# %%
