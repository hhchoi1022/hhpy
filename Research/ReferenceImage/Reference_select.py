
"""
Spyder Editor

This is a temporary script file.
"""
from astropy.io import fits
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import pandas as pd

#%%

observatory = 'KCT_STX16803'
folder = f'/data3/hhchoi1022/IMSNG/{observatory}'
targetlist= sorted(os.listdir(f'{folder}/Data'))
if observatory == 'KCT_STX16803':
    filterlist =  {'g':'green','r':'red','i':'black'}
if observatory == 'CBNUO':
    filterlist = {'B':'blue','V':'green','R':'red','I':'black'}



for filter_ in filterlist:
    depth_R=[]
    seeing_R=[]
    color = filterlist[filter_]
    for target in targetlist:
        imkey = f'{folder}/Data/{target}/{observatory}/{filter_}/*.fits'
        imagelist = glob.glob(imkey)
        for image in imagelist:
            hdr = fits.getheader(image)
            if 'HSEEING' in hdr.keys():
                seeing_R.append(hdr['HSEEING'])
                depth_R.append(hdr['HDEPTH5'])
    cut_seeing_R_1 = round(np.median(seeing_R)-0.7*np.std(seeing_R),1)
    cut_seeing_R_2 = round(np.median(seeing_R),1)
    cut_seeing_R_3 = round(np.median(seeing_R)+0.7*np.std(seeing_R),1)
    cut_seeing_R_4 = round(np.median(seeing_R)+3*np.std(seeing_R),1)
    cut_depth_R_1 = round(np.median(depth_R)+0.5*np.std(depth_R),3)
    cut_depth_R_2 = round(np.median(depth_R),3)
    cut_depth_R_3 = round(np.median(depth_R)-0.5*np.std(depth_R),3)
    cut_depth_R_4 = round(np.median(depth_R)-3*np.std(depth_R),3)
    
    plt.xlabel('Seeing[arcsec]')
    plt.ylabel('Depth[AB]')
    plt.plot(seeing_R, depth_R, marker='o', mec=color, mfc='none', ls='none', label=f'{filter_} Filter[{len(seeing_R)}]', alpha=0.1)
    plt.axvline(x=cut_seeing_R_1, c='green', ls='-', label=f'Seeing <= {cut_seeing_R_1}')
    plt.axhline(y=cut_depth_R_1, c='green', ls='-', label=f'Depth >= {cut_depth_R_1}')
    plt.axvline(x=cut_seeing_R_2, c='black', ls='-', label=f'Seeing <= {cut_seeing_R_2}')
    plt.axhline(y=cut_depth_R_2, c='black', ls='-', label=f'Depth >= {cut_depth_R_2}')
    plt.axvline(x=cut_seeing_R_3, c='blue', ls='-', label=f'Seeing <= {cut_seeing_R_3}')
    plt.axhline(y=cut_depth_R_3, c='blue', ls='-', label=f'Depth >= {cut_depth_R_3}')
    plt.axvline(x=cut_seeing_R_4, c='red', ls='-', label=f'Seeing <= {cut_seeing_R_4}')
    plt.axhline(y=cut_depth_R_4, c='red', ls='-', label=f'Depth >= {cut_depth_R_4}')
    plt.legend()
    os.makedirs(f'{folder}/selected',exist_ok = True)
    os.chdir(f'{folder}/selected')
    plt.savefig(f'{observatory}_{filter_}_cutline.png') ######
    plt.show()




    best_R = []
    normal_R = []
    bad_R = []
    worst_R = []
    summary = pd.DataFrame(index = targetlist, columns = ['#tot','#sel','status'])

    for target in targetlist:
        depth_R=[]
        seeing_R=[]
        best = []
        normal = []
        bad =[]
        worst =[]
        imkey = f'{folder}/Data/{target}/{observatory}/{filter_}/*.fits'
        imagelist = glob.glob(imkey)
        color = filterlist[filter_]
        n_tot = len(imagelist)
        for image in imagelist:
            hdr = fits.getheader(image)
            if 'HSEEING' in hdr.keys():
                seeing_R.append(hdr['HSEEING'])
                depth_R.append(hdr['HDEPTH5'])
                if hdr['HSEEING'] <cut_seeing_R_1 and hdr['HDEPTH5']>cut_depth_R_1:
                    best.append(image)
                    best_R.append(image)
                if hdr['HSEEING'] <cut_seeing_R_2 and hdr['HDEPTH5'] >cut_depth_R_2:
                    normal_R.append(image)
                    normal.append(image)
                if hdr['HSEEING'] < cut_seeing_R_3 and hdr['HDEPTH5'] > cut_depth_R_3:
                    bad_R.append(image)
                    bad.append(image)
                if hdr['HSEEING'] < cut_seeing_R_4 and hdr['HDEPTH5'] > cut_depth_R_4:
                    worst_R.append(image)
                    worst.append(image)
            
        
        path = f'/{folder}/selected/{target}/{filter_}'
        os.makedirs(path,exist_ok=True)
        os.chdir(path)
        if len(best)>=20:
            status = 'best'
            n_sel = len(best)
            os.makedirs('best',exist_ok=True)
            for image in best: 
                os.system(f'cp {image} {path}')
        elif len(normal)>=20:
            status = 'normal'
            n_sel = len(normal)
            os.makedirs('normal',exist_ok=True)
            for image in normal:
                os.system(f'cp {image} {path}')
        elif len(bad)>=10:
            status = 'bad'
            n_sel = len(bad)
            os.makedirs('bad',exist_ok=True)
            for image in bad:
                os.system(f'cp {image} {path}')   
        else:
            status = 'worst'
            n_sel = len(worst)
            os.makedirs('worst',exist_ok=True)
            for image in worst:
                 os.system(f'cp {image} {path}')
        
        summary['status'][target] = status
        summary['#tot'][target] = n_tot
        summary['#sel'][target] = n_sel
        plt.figure(figsize = (7,5))
        plt.title(f'{target}/{filter_}')
        plt.xlabel('Seeing[arcsec]')
        plt.ylabel('Depth[AB]')
        plt.xlim(cut_seeing_R_1-1,cut_seeing_R_4+1)
        plt.ylim(cut_depth_R_4-1,cut_depth_R_1+1)
        plt.plot(seeing_R, depth_R, marker='o', mec=color, mfc='none', ls='none', label=f'{filter_} Filter[{len(seeing_R)}]', alpha=0.3)
        plt.axvline(x=cut_seeing_R_1, c='green', ls='-',alpha = 0.5)
        plt.axhline(y=cut_depth_R_1, c='green', ls='-', label=f'Best[{len(best_R)}]',alpha = 0.5)
        plt.axvline(x=cut_seeing_R_2, c='black', ls='-',alpha = 0.5)
        plt.axhline(y=cut_depth_R_2, c='black', ls='-', label=f'Normal[{len(normal_R)}]', alpha = 0.5)
        plt.axvline(x=cut_seeing_R_3, c='red', ls='-',alpha = 0.5)
        plt.axhline(y=cut_depth_R_3, c='red', ls='-', label=f'Bad[{len(bad_R)}]', alpha = 0.5)
        plt.axvline(x=cut_seeing_R_4, c='red', ls='-',alpha = 0.5)
        plt.axhline(y=cut_depth_R_4, c='red', ls='-', label=f'Worst[{len(worst_R)}]', alpha = 0.5)
        plt.legend(loc = 3)
        os.chdir(f'/{folder}/selected/{target}')
        plt.savefig(f'{target}_{filter_}.png')
        plt.show()
    os.chdir(f'{folder}/selected')
    summary.to_csv(f'summary_{filter_}.csv')
#%%
