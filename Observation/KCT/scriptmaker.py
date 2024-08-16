#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 12:35:09 2021

@author: Mankeun Jeong
Last Update: 2021.05.04.

"""
#%%
import time
import copy
import calendar
import os, glob
import numpy as np
import pandas as pd
import datetime as dt
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table
import rtsmaker
from astropy.coordinates import SkyCoord
#%%
def scriptmaker(obs, year, month, day, exp=60, cnt=12, temp=-15.0, path_input='./', path_output='./'):

    #0. time parameter
    
    if obs=='RASA36' or 'KCT':
        LT2UTC      = dt.timedelta(hours=3) # Chilian local time to UTC
    else:
        LT2UTC      = dt.timedelta(hours=int(input("Enter the time difference btw LT and UTC:")))
    TOMORROW    = dt.timedelta(days=1)
    
    #1. Reading tonight's RTS file from the input directory.
    
    project     = 'IMSNG'
    
    os.chdir(path_input)
    inputfile = '{}-{}{}{}-rts_vis-{}.txt'.format(project, year, month, day, obs)
    rts         = ascii.read(inputfile)
    rts.sort('ra')

    #2. Checking the twilight time. (Starting & ending of an observation)
    
    with open(inputfile, 'r') as f:
        contents    = f.readlines()
        for c in contents:
            if 'deg sunset' in c:
                sunset  = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, c.split()[-1]), '%Y-%m-%d %H:%M') + LT2UTC
            elif 'deg sunrise' in c:
                sunrise = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, c.split()[-1]), '%Y-%m-%d %H:%M') + TOMORROW + LT2UTC
                
    #3. Writing the script.
    
    os.chdir(path_output)
    with open('script_{}{}{}_{}_{}.txt'.format(year, month, day, project, obs), 'w') as f:
        
        # starting time
        f.write('; Prepare for the observation\n')
        prepare         = sunset-dt.timedelta(minutes=30)
        prepare_form    = '{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(prepare.month, prepare.day, prepare.year%2000, prepare.hour, prepare.minute, prepare.second)
        sunset_form     = '{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(sunset.month, sunset.day, sunset.year%2000, sunset.hour, sunset.minute, sunset.second)
        sunrise_form    = '{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(sunrise.month, sunrise.day, sunrise.year%2000, sunrise.hour, sunrise.minute, sunrise.second)
        f.write('#WAITUNTIL 1, {}\n'.format(prepare_form))
    
        # cooling
        f.write('\n; Cooling\n')
        f.write('#CHILL {}, {}\n'.format(temp, 1.0))
        
        # calibration
        f.write('\n; Calibration frames\n')
        f.write("#COUNT {}\n#INTERVAL {}\n#DARK\n".format(cnt, exp)) # dark
        f.write("#COUNT {}\n#BIAS\n".format(cnt)) # bias
        
        f.write('\n; Start of the evening\n')
        f.write('#WAITUNTIL 1, {}\n'.format(sunset_form))

        # autofocus
        f.write('\n; Autofocus\n')
        f.write('#AFINTERVAL {}\n'.format(180)) # AF every 3 hrs
    
        # targeting
        f.write('\n; Targeting')
    
        total       = 0
        midnight    = 0
        
        for i in range(len(rts)):

            if i !=0:
                dumtime1    = int(rts['transit(LT)'][i-1][:2])
                dumtime2    = int(rts['transit(LT)'][i][:2])
                if dumtime1 - 20 > dumtime2: # date passed (ex. 23:59 --> 00:01)
                    midnight += 1
     
            target  = rts['name'][i]
            ra      = rts['ra'][i]
            dec     = rts['dec'][i]
            
            if midnight == 0:
                transit     = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, rts['transit(LT)'][i]), '%Y-%m-%d %H:%M') + LT2UTC
                rise        = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, rts['rise(LT)'][i]), '%Y-%m-%d %H:%M') + LT2UTC
            elif midnight == 1:
                transit     = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, rts['transit(LT)'][i]), '%Y-%m-%d %H:%M') + TOMORROW + LT2UTC
                if int(rts['rise(LT)'][i][:2]) < sunrise.hour:
                    rise    = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, rts['rise(LT)'][i]), '%Y-%m-%d %H:%M') + TOMORROW + LT2UTC
                else: #(ex. rises at 22:00 & transits at 01:00)
                    rise    = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, rts['rise(LT)'][i]), '%Y-%m-%d %H:%M') + LT2UTC
            else:
                print('Date passed twice. Check the file or parameters.')
            
            if transit - rise > dt.timedelta(days=1):
                transit     = transit - dt.timedelta(days=1)
            
            # WAITTIME = Average of TRANSITTIME and RISINGTIME. (22.01.18. update)
            # For more description please check the "WAITUNTIL" directive from: http://solo.dc3.com/ar/RefDocs/HelpFiles/ACP81Help/planfmt.html#directives
            
            wait       = transit - (transit - rise)/2
            wait_form  = '{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(wait.month, wait.day, wait.year%2000, wait.hour, wait.minute, wait.second)
            """
            if (transit - dt.timedelta(hours=2)) > rise:
                wait       = transit - dt.timedelta(hours=2)
                if wait > sunrise:
                    wait    = rise
                wait_form  = '{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(wait.month, wait.day, wait.year%2000, wait.hour, wait.minute, wait.second)
            else:
                wait       = rise
                wait_form  = '{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(wait.month, wait.day, wait.year%2000, wait.hour, wait.minute, wait.second)
            """
            
            if wait < sunrise:
                f.write('\n#WAITUNTIL 1, {}\n'.format(wait_form))
                f.write('#COUNT {}\n#INTERVAL {}\n#NOPREVIEW\n#NOSOLVE\n{}\t{}\t{}\n'.format(cnt, exp, target, ra, dec))
                total   = total + 1
            elif rise < sunrise:
                wait = rise
                wait_form  = '{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(wait.month, wait.day, wait.year%2000, wait.hour, wait.minute, wait.second)
                f.write('\n#WAITUNTIL 1, {}\n'.format(wait_form))
                f.write('#COUNT {}\n#INTERVAL {}\n#NOPREVIEW\n#NOSOLVE\n{}\t{}\t{}\n'.format(cnt, exp, target, ra, dec))
                total   = total + 1
            else:
                print("{} rises after moring".format(target))
                pass # Targets rising after the sunrise will not be included.
    
        # closing

        f.write('#QUITAT {}\n'.format(sunrise_form))

    return 0

def sorting_scripts(path_script):
    
    os.chdir(path_script)
    allscript = sorted(glob.glob('script*.txt'))
    
    for i in range(len(allscript)):
        # target sorting
        f = open(allscript[i], 'r').readlines()
    
        seg     = []
        for j in range(len(f)):
            if f[j].startswith('#WAITUNTIL 1') and f[j+1].startswith('#COUNT'):
                    seg.append('{}{}{}{}{}{}'.format(f[j], f[j+1], f[j+2], f[j+3], f[j+4], f[j+5]))
        
        newseg = sorted(seg)
        # rewrite
        with open('{}'.format(allscript[i]), 'w') as out:
            
            add = True
            for line in f:
                if line.startswith('#AFINTERVAL'):
                    out.write(line)
                    out.write('\n;Targeting\n')
                    add     = False
                else:
                    if add:
                        out.write(line)
                    else:
                        pass
                    
            for newline in newseg:
                out.write(newline+'\n')
            
            add = False
            for line in f:
                if line.startswith('; Closing'):
                    out.write(';Total Targets: {}\n\n'.format(len(seg)))
                    out.write(line)
                    add     = True
                else:
                    if add:
                        out.write(line)
                    else:
                        pass
    return 0

def add_new_target(name, ra, dec, y, m, d, path_tbl, path_script, obs='RASA36'):
    """
    name   = 'SDSSJ1430+2303'
    ra     = '14:30:16.05' 
    dec    = '+23:03:44.4'
    y,m,d  = 2022,3,6
    """
    # basic information
    observer    = rtsmaker.callobserver(obs, ascii.read(f'{path_tbl}/observatory.txt'))
    obsdat      = Table([[name], [ra], [dec], [1]], names=['obj', 'ra', 'dec', 'priority'])
    coord       = SkyCoord(obsdat['ra'], obsdat['dec'], unit=(u.hourangle, u.deg))
    
    # rts generator
    rtsmaker.rts_vis_maker(coord, observer, y, m, d, obs=obs, intbl=obsdat, constraint_alt=30*u.deg, constraint_moon=40*u.deg, path_out=path_tbl)
    
    # script generator
    scriptmaker(obs, str(y).zfill(4), str(m).zfill(2), str(d).zfill(2), exp=60, cnt=12, temp=-13.0, path_input=path_tbl, path_output=path_tbl)
    
    # update the original script
    f_origin    = open(f'{path_script}/script_{y:04d}{m:02d}{d:02d}_IMSNG_{obs}.txt', 'r').readlines()
    f_update    = open(f'{path_tbl}/script_{y:04d}{m:02d}{d:02d}_IMSNG_{obs}.txt', 'r').readlines()
    
    seg     = []
    for j in range(len(f_update)):
        if f_update[j].startswith('#WAITUNTIL 1') and f_update[j+1].startswith('#COUNT'):
                seg.append('{}{}{}{}{}{}'.format(f_update[j], f_update[j+1], f_update[j+2], f_update[j+3], f_update[j+4], f_update[j+5]))
    
    with open(f'{path_script}/script_{y:04d}{m:02d}{d:02d}_IMSNG_{obs}_updated.txt', 'w') as out:
        
        for line in f_origin:
            if line.startswith(';Total Target'):
                out.write(seg[0]+'\n') # new target script
                out.write(line+'\n')
                break
            else:
                out.write(line)
        out.write(f_origin[-2]) # ; Closing
        out.write(f_origin[-1]) # it should be sth like '#QUITAT 03/07/22 09:15:00\n'
    
    # obs time optimization
    sorting_scripts(path_script)

    return 0

#%% USER SETTING
obs         = 'KCT'
project     = 'IMSNG'

# pathes
path_base   = '/home/sonic/research/observation/script'
path_obs    = '/home/sonic/research/observation/script/observatory.txt'
path_rts    = '/home/sonic/research/observation/script/rts_vis_{}'.format(2022)
path_script = '/home/sonic/research/observation/script/script_{}'.format(2022)
path_tbl    = '/home/sonic/research/observation/script/target_table'
# path_mscript= '/home/sonic/research/observation/script/script_{}/merged'.format(2022)
path_rts = '/home/hhchoi1022/Desktop/Gitrepo/observatory/KCT'
path_script = './'
#%% clearance
clean_script    = True

if clean_script:
    os.chdir(path_script)
    scripts     = sorted(glob.glob('*.txt'))
    for s in scripts:
        os.system('rm {}'.format(s))

# clean_mscript   = False
# if clean_mscript:
#     os.chdir(path_mscript)
#     mscripts     = sorted(glob.glob('*.txt'))
#     for m in mscripts:
#         os.system('rm {}'.format(m))

os.chdir(path_base)    
#%% Daily script maker
now         = dt.datetime.now() - dt.timedelta(hours=12) # for Chilean Local Time
year        = str(now.year)
month       = str(now.month).zfill(2)
day         = str(now.day).zfill(2)

scriptmaker(obs, year, month, day, exp=60, cnt=12, temp=-13.0, path_input=path_rts, path_output=path_script)
#%% Monthly script maker
for date in calendar.Calendar().itermonthdates(2022,4):

    year    = str(date.year)
    month   = str(date.month).zfill(2)
    day     = str(date.day).zfill(2)

    scriptmaker(obs, year, month, day, exp=60, cnt=12, temp=-13.0, path_input=path_rts, path_output=path_script)
#%% sorting mscripts

os.chdir(path_script)
allscript = sorted(glob.glob('script*.txt'))

for i in range(len(allscript)):
    # target sorting
    f = open(allscript[i], 'r').readlines()

    seg     = []
    for j in range(len(f)):
        if f[j].startswith('#WAITUNTIL 1') and f[j+1].startswith('#COUNT'):
                seg.append('{}{}{}{}{}{}'.format(f[j], f[j+1], f[j+2], f[j+3], f[j+4], f[j+5]))
    
    newseg = sorted(seg)
    # rewrite
    with open('{}'.format(allscript[i]), 'w') as out:
        
        add = True
        for line in f:
            if line.startswith('#AFINTERVAL'):
                out.write(line)
                out.write('\n;Targeting\n')
                add     = False
            else:
                if add:
                    out.write(line)
                else:
                    pass
                
        for newline in newseg:
            out.write(newline+'\n')
        
        add = False
        for line in f:
            if line.startswith('; Closing'):
                out.write(';Total Targets: {}\n\n'.format(len(seg)))
                out.write(line)
                add     = True
            else:
                if add:
                    out.write(line)
                else:
                    pass
    

#%% Script merger: No point
# os.chdir(path_script)
# allscript   = glob.glob('script*.txt'); allscript.sort()

# for i in range(len(allscript)-1):
    
#     with open('merged/m{}'.format(allscript[i]), 'w') as out:
        
#         tonight     = open(allscript[i], 'r').readlines()
#         tomorrow    = open(allscript[i+1], 'r').readlines()
        
#         for line1 in tonight:
#             if line1.startswith('; Closing'):
#                 add    = False
#                 for line2 in tomorrow:
#                     if line2.startswith('; Targeting'):
#                         out.write(line2)
#                         add     = True
#                     elif line2.startswith('; Closing'):
#                         break
#                     else:
#                         if not(add):
#                             pass
#                         else:
#                             out.write(line2)
#                 out.write(line1)
#             else:
#                 out.write(line1)
        