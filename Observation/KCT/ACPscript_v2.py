#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 17:19:36 2022

@author: hhchoi1022
"""
"""
conf_param = dict(
    #obspath         = "/home/hhchoi1022/Desktop/Gitrepo/Config/obs_location.txt",
    #catpath         = '/home/hhchoi1022/Desktop/Gitrepo/Config/alltarget.dat'
    rtspath = '/home/hhchoi1022/Desktop/Gitrepo/observatory/KCT/ACP/IMSNG-20170311-rts_vis-KCT.txt',
    sptpath = '/home/hhchoi1022/Desktop/Gitrepo/observatory/KCT/',
    obs = 'KCT',
    project = 'IMSNG',x
    grid_time = 8000,
    exptime = 1200,
    ntarget_grid = 6#, # Number of targets in each grid
    #starttime = '2022/04/22 01:38:00'
    )
"""
#%%

import datetime
import datetime as dt
import os
import re

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pytz
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.io import ascii
from astropy.table import Table, vstack
from astropy.time import Time
from HHsupport_phot import to_skycoord

from rtsmaker import rtsmaker


class ACP_scriptmaker:
    
    def __init__(self, 
                 rtspath,
                 calib_interval = 10800,
                 overhead_ratio = 1.10,
                 select_mode = 'optimal',
                 observatory = 'KCT',
                 obsinfopath = "/home/hhchoi1022/Desktop/Gitrepo/Config/obs_location.txt",
                 savefigpath = None,
                 scriptpath = None,
                 safe_time = 1200,
                 weight_altitude = 1,
                 weight_priority = 1,
                 ccdtemp = -10,
                 dark = True,
                 bias = True,
                 flat = True
                 ):
        ####################################
        ###### Setting observatory
        ####################################
        obsinfo               = ascii.read(obsinfopath)
        obsname               = np.copy(obsinfo['name'])
        obsindex              = np.where(obsname == observatory)[0]
        obslat                = (np.copy(obsinfo['latitude(N+)'])[obsindex])[0]
        obslon                = (np.copy(obsinfo['longitude(E+)'])[obsindex])[0]
        obsalt                = (np.copy(obsinfo['altitude'])[obsindex])[0]
        obstz                 = (np.copy(obsinfo['timezone'])[obsindex])[0]
        self.rtspath          = rtspath
        self.tz               = pytz.timezone(obstz)
        self.loc              = EarthLocation(lat=obslat*u.deg, lon=obslon*u.deg, height=obsalt*u.m)
        self.overhead_ratio   = overhead_ratio
        self.select_mode      = select_mode
        self.safe_time        = safe_time
        self.savefigpath      = savefigpath
        self.calib_interval   = calib_interval
        self.weight_altitude  = weight_altitude
        self.weight_priority  = weight_priority
        self.rts, self.sunset, self.sunrise, self.LT2UTC = self.read_rts(rtspath = rtspath)
        
        self.timegrid         = self.set_timegrid()

        self.all_targetlist   = self.construct_plan()
        if scriptpath == None:
            scriptpath = os.path.dirname(rtspath)+'/ACP_'+os.path.basename(rtspath)
        if savefigpath == None:
            self.savefigpath = os.path.dirname(rtspath)+'/ACP_'+os.path.basename(rtspath)[:-3]+'png'
        self.show_graph()
        self.scriptpath       = self.make_ACPscript(filename = scriptpath, ccdtemp = ccdtemp, dark = dark, bias = bias, flat = flat)

    def read_rts(self, rtspath):
        ####################################
        ###### Reading rts information 
        ####################################
        rts         = ascii.read(rtspath)
        rts['name_tmp'] = '                '
        rts['exptime_tmp'] = '                '
        rts['filter_tmp'] = '                '
        rts['counts_tmp'] = '                '
        rts['binning_tmp'] = '                '
        rts.sort('rise(LT)')
        with open(rtspath, 'r') as f:
            contents    = f.readlines()
            year, month, day = re.findall('(\d{4})/(\d{1,2})/(\d{1,2})',contents[0])[0]
            year = int(year)
            month = int(month)
            day = int(day) 
            datetime_now = datetime.datetime(year,month,day)-datetime.timedelta(days = 1)
            year = datetime_now.year
            month = datetime_now.month
            day = datetime_now.day
            LT2UTC = dt.timedelta(hours = -self.tz.utcoffset(datetime_now, is_dst = True).total_seconds()/3600)
            TOMORROW    = dt.timedelta(days=1)
            for c in contents:
                if 'deg sunset' in c:
                    sunsettime  = re.findall(r'(\d{1,2}:\d{1,2})',c)[0]
                    sunset = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, sunsettime), '%Y-%m-%d %H:%M') + LT2UTC
                elif 'deg sunrise' in c:
                    sunrisetime  = re.findall(r'(\d{1,2}:\d{1,2})',c)[0]
                    sunrise = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, sunrisetime), '%Y-%m-%d %H:%M') + TOMORROW + LT2UTC
        rdtimelist = []
        tdtimelist = []
        sdtimelist = []
        for i in range(len(rts)):
            rtime = [int(re.findall('(\d{2}):(\d{2})',rts['rise(LT)'][i])[0][0]), int(re.findall('(\d{2}):(\d{2})',rts['rise(LT)'][i])[0][1])]
            ttime = [int(re.findall('(\d{2}):(\d{2})',rts['transit(LT)'][i])[0][0]), int(re.findall('(\d{2}):(\d{2})',rts['transit(LT)'][i])[0][1])]
            stime = [int(re.findall('(\d{2}):(\d{2})',rts['set(LT)'][i])[0][0]), int(re.findall('(\d{2}):(\d{2})',rts['set(LT)'][i])[0][1])]
        
            # date passed 23:59 >> 00:01 
            midnight = 0
            if (i < int(len(rts)/2)) & (rtime[0] < 10):
                midnight = 1
            
            if midnight == 0 :
                rdtime = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, rts['rise(LT)'][i]), '%Y-%m-%d %H:%M') + LT2UTC
                tdtime = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, rts['transit(LT)'][i]), '%Y-%m-%d %H:%M') + LT2UTC
                sdtime = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, rts['set(LT)'][i]), '%Y-%m-%d %H:%M') + LT2UTC
            else:
                rdtime = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, rts['rise(LT)'][i]), '%Y-%m-%d %H:%M') + LT2UTC + TOMORROW
                tdtime = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, rts['transit(LT)'][i]), '%Y-%m-%d %H:%M') + LT2UTC + TOMORROW
                sdtime = dt.datetime.strptime('{}-{}-{} {}'.format(year, month, day, rts['set(LT)'][i]), '%Y-%m-%d %H:%M') + LT2UTC + TOMORROW
            
            if (sdtime - rdtime).days <0:
                sdtime += TOMORROW
            if (tdtime - rdtime).days <0:
                tdtime += TOMORROW
            
            rdtimelist.append(rdtime)
            tdtimelist.append(tdtime)
            sdtimelist.append(sdtime)
        rts['name_tmp'][:] = rts['name']
        rts['exptime_tmp'][:] = rts['exptime']
        rts['filter_tmp'][:] = rts['filter']
        rts['counts_tmp'][:] = rts['counts']
        rts['binning_tmp'][:] = rts['binning']
        rts.remove_columns(['name', 'exptime', 'filter','counts','binning'])
        rts.rename_columns(['name_tmp', 'exptime_tmp', 'filter_tmp','counts_tmp','binning_tmp'],['name', 'exptime', 'filter','counts','binning'])
        rts['rise(UT)'] = rdtimelist
        rts['transit(UT)'] = tdtimelist
        rts['set(UT)'] = sdtimelist
        return rts, sunset, sunrise, LT2UTC 

    def set_timegrid(self):
        timegrid = []
        time = self.sunset
        while time < self.sunrise:
            timegrid.append(time)
            time = time + dt.timedelta(seconds = self.calib_interval) 
        timegrid.append(self.sunrise)
        return timegrid

    def select_target(self,
                      rts_tbl,
                      time,
                      mode = 'altitude',
                      safe_time = 1800,
                      weight_altitude = 1,
                      weight_priority = 1
                      ):
        tmp_tbl = rts_tbl.copy()
        alt_score = weight_altitude * ((np.abs(tmp_tbl['transit(UT)'] - time)) / (tmp_tbl['set(UT)'] - tmp_tbl['rise(UT)']))
        #timeangle = (tmp_tbl['transit(UT)'] - time)
        #timeangle[(timeangle < 0) & (timeangle > -safe_time)]
        alt_score = weight_altitude * ((tmp_tbl['transit(UT)'] - time) / (tmp_tbl['set(UT)'] - tmp_tbl['rise(UT)']))
        priority_score = weight_priority * tmp_tbl['priority']/np.max(tmp_tbl['priority'])
        tmp_tbl['alt_score'] = alt_score
        tmp_tbl['priority_score'] = priority_score
        tmp_tbl['tot_score'] = np.mean([alt_score, priority_score], axis=  0)
        tmp_tbl['Until_set'] = tmp_tbl['set(UT)'] - time
        tmp_tbl = tmp_tbl[(tmp_tbl['rise(UT)'] < (time - dt.timedelta(seconds =safe_time))) & (tmp_tbl['set(UT)'] > (time + dt.timedelta(seconds =safe_time)))]

        if mode.lower() =='urgent':
            if len(tmp_tbl) > 0:
                tmp_tbl['exp_start(UT)'] = time
                tmp_tbl.sort('Until_set')
                return tmp_tbl[0]
            else:
                return Table()
            
        elif mode.lower() == 'optimal':
            if len(tmp_tbl) > 0:
                tmp_tbl['exp_start(UT)'] = time
                tmp_tbl.sort('tot_score')
                return tmp_tbl[0]
            else:
                return Table()
            
        elif mode.lower() == 'altitude':
            if len(tmp_tbl) > 0:
                tmp_tbl['exp_start(UT)'] = time
                tmp_tbl.sort(['alt_score','priority_score'])
                return tmp_tbl[0]
            else:
                return Table()
            
        elif mode.lower() == 'priority':
            if len(tmp_tbl) > 0:
                tmp_tbl['exp_start(UT)'] = time
                tmp_tbl.sort(['priority_score','alt_score'])
                return tmp_tbl[0]
            else:
                return Table()
        
    def calc_totexp(self,
                    target):
        tot_exptime = np.sum([int(exp)*int(count) for exp, count in zip(target['exptime'].split(','), target['counts'].split(','))])
        return tot_exptime
    
    def construct_plan(self):
        # for the first time slot, select the urgent targets before the target set
        all_targetlist = Table()
        targetlist = Table()
        obstime = self.timegrid[0]
        remain_tbl = self.rts
        while obstime + dt.timedelta(seconds = self.safe_time) < self.timegrid[1]:
            best_target = self.select_target(remain_tbl, obstime, mode = 'urgent')
            if len(best_target) >0:
                remain_tbl = remain_tbl[remain_tbl['name'] != best_target['name']]
                tot_exptime = self.calc_totexp(best_target)
                obstime += dt.timedelta(seconds = tot_exptime * self.overhead_ratio)
            else:
                remain_tbl = remain_tbl
                tot_exptime = self.calc_totexp(remain_tbl[0])
                obstime += dt.timedelta(seconds = tot_exptime * self.overhead_ratio)
            targetlist = vstack([best_target, targetlist])
        targetlist.sort('exp_start(UT)')
        targetlist['group'] = 0 
        all_targetlist = vstack([targetlist, all_targetlist])
        for i in range(len(self.timegrid)):
            if (i != 0) & (i != len(self.timegrid)-1):
                obstime = self.timegrid[i] 
                next_timegrid = self.timegrid[i+1]
                targetlist = Table()
                while (obstime + dt.timedelta(seconds = self.safe_time) < next_timegrid) & (len(remain_tbl) > 0):
                        best_target = self.select_target(remain_tbl, obstime, mode = self.select_mode)
                        if len(best_target) >0:
                            remain_tbl = remain_tbl[remain_tbl['name'] != best_target['name']]
                            tot_exptime = self.calc_totexp(best_target)
                            obstime += dt.timedelta(seconds = tot_exptime * self.overhead_ratio )
                        else:
                            remain_tbl = remain_tbl
                            tot_exptime = self.calc_totexp(remain_tbl[0])
                            obstime += dt.timedelta(seconds = tot_exptime * self.overhead_ratio)
                        targetlist = vstack([best_target, targetlist])
                if len(targetlist) > 0:
                    targetlist.sort('exp_start(UT)')
                    targetlist['group'] = i
                    all_targetlist = vstack([targetlist, all_targetlist])
        all_targetlist.sort('exp_start(UT)')

        # Priority 0 target must be observed just after rising
        prior0list = self.rts[self.rts['priority'] == 0.0]
        for prior0 in prior0list:
            if (prior0['name'] not in all_targetlist['name']):
                closest_idx = np.argmin(np.abs(all_targetlist['exp_start(UT)']-(prior0['rise(UT)']+dt.timedelta(seconds = 1800))))
                for colname in prior0.colnames:
                    all_targetlist[closest_idx][colname] = prior0[colname]
                starttime = all_targetlist[closest_idx]['exp_start(UT)']
                tot_exptime = self.calc_totexp(all_targetlist[closest_idx])
                endtime = starttime + dt.timedelta(seconds = int(tot_exptime))
                remove_idx = []
                for i, startend in enumerate(zip(all_targetlist['exp_start(UT)']-starttime, all_targetlist['exp_start(UT)']-endtime)):
                    start = startend[0]
                    end = startend[1]
                    if (start.total_seconds() > 0) & (end.total_seconds() <= 0):
                        remove_idx.append(i)
                print('--'*60)
                all_targetlist.remove_rows(remove_idx)
            if (prior0['name'] == all_targetlist['name'][-1]):
                closest_idx = np.argmin(np.abs(all_targetlist['exp_start(UT)']-(prior0['rise(UT)']+dt.timedelta(seconds = 1800))))
                for colname in prior0.colnames:
                    all_targetlist[closest_idx][colname] = prior0[colname]
                starttime = all_targetlist[closest_idx]['exp_start(UT)']
                tot_exptime = self.calc_totexp(all_targetlist[closest_idx])
                endtime = starttime + dt.timedelta(seconds = int(tot_exptime))
                remove_idx = []
                for i, startend in enumerate(zip(all_targetlist['exp_start(UT)']-starttime, all_targetlist['exp_start(UT)']-endtime)):
                    start = startend[0]
                    end = startend[1]
                    if (start.total_seconds() > 0) & (end.total_seconds() < 0):
                        remove_idx.append(i)
                all_targetlist.remove_rows(remove_idx)
                best_target = self.select_target(remain_tbl, all_targetlist['exp_start(UT)'][-1], mode = self.select_mode)
                for colname in best_target.colnames:
                    all_targetlist[-1][colname] = best_target[colname]
            print(60*'=')
        return all_targetlist

    def show_graph(self):
        t = Time(self.sunset,scale='utc')
        LT2UTC = dt.timedelta(hours = -self.tz.utcoffset(t.datetime, is_dst = True).total_seconds()/3600)
        delta_midnight = np.linspace(-2, 12, 300)*u.hour
        delta_timegrid = [(self.timegrid[i]-self.timegrid[0]).seconds/3600 for i in range(len(self.timegrid))]
        times = t+delta_midnight
        earthframe = AltAz(obstime=times, location=self.loc)

        plt.figure(figsize = (15,6), dpi = 300)
        plt.rcParams["figure.figsize"] = (15,6)
        plt.ylim(10,90)
        plt.xlim(-2, 12)
        plt.fill_betweenx([-90,90], 0, (self.sunrise-self.sunset).seconds/3600, color = 'k', alpha = 0.1)
        plt.xlabel(f'UT - sunset [{self.sunset-LT2UTC}]', fontsize=20)
        plt.ylabel('Altitude [deg]', fontsize=20)
        plt.axhline(30, linestyle = '--', c = 'r', label = 'altitude limit')
        for timeslot in delta_timegrid:
            plt.axvline(timeslot, linestyle = '--', c = 'k', alpha = 0.3)
        plt.text(-1.8, 75, f'Mode = {self.select_mode}')
        plt.text(-1.8, 72, f'# of target = {len(self.all_targetlist)}')
        for target in self.all_targetlist:
            ra, dec = to_skycoord(target['ra'],target['dec']).ra.value, to_skycoord(target['ra'],target['dec']).dec.value
            coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5')
            altaz = coord.transform_to(earthframe)
            
            obsstart_time = target['exp_start(UT)']
            obsend_time = target['exp_start(UT)'] + dt.timedelta(seconds = self.calc_totexp(target) * self.overhead_ratio )
            obs_idx = np.where([(times >= obsstart_time) & (times <= obsend_time)])[1]
            obs_times = delta_midnight[obs_idx]
            obs_altaz = altaz[obs_idx]
            plt.plot(delta_midnight, altaz.alt, c = 'k', alpha = 0.2)
            plt.plot(obs_times, obs_altaz.alt, c = 'k')
            if target['priority'] == 0.0:
                plt.plot(obs_times, obs_altaz.alt, c = 'b')
            if target['priority'] == 5.0:
                ToO_ra, ToO_dec = to_skycoord(target['ra'],target['dec']).ra.value, to_skycoord(target['ra'],target['dec']).dec.value
                ToO_coord = SkyCoord(ra=ToO_ra*u.degree, dec=ToO_dec*u.degree, frame='fk5')
                ToO_altaz = ToO_coord.transform_to(earthframe)
                obs_ToO_times = delta_midnight[obs_idx]
                obs_ToO_altaz = ToO_altaz[obs_idx]
                plt.axvline((target['exp_start(UT)']-self.timegrid[0]).seconds/3600, linestyle = '--', c = 'r', alpha = 1)
                plt.annotate('ToO start', ha = 'center', va = 'bottom',xytext = (np.abs((target['exp_start(UT)']-self.timegrid[0])).seconds/3600-1, 85),
                xy = (np.abs((target['exp_start(UT)']-self.timegrid[0])).seconds/3600+2, 85), c= 'r',
                arrowprops = { 'facecolor' : 'r', 
                            'edgecolor':'k', 
                            'shrink' : 0.2, 
                            'alpha':1
                            })
                plt.plot(delta_midnight, ToO_altaz.alt, c = 'r', alpha = 0.2)
                plt.plot(obs_ToO_times, obs_ToO_altaz.alt, c = 'r')
        plt.plot(-99, -99, c = 'b', label = 'prioirty 0')     
        plt.plot(-99, -99, c = 'k', label = 'observation')
        plt.plot(-99, -99, c = 'r', label = 'ToO')
        plt.legend(loc = 1)
        if self.savefigpath != None:
            plt.savefig(self.savefigpath)
    
    def make_ACPscript(self,
                      filename = './plan.txt',
                      ccdtemp = -10,
                      dark = True,
                      bias = True,
                      flat = True
                      ):
        
        def dt_to_form(dt):
            form ='{:02}/{:02}/{:02} {:02}:{:02}:{:02}'.format(dt.month, dt.day, dt.year%2000, dt.hour, dt.minute, dt.second)
            return form
        def dt_to_form_flat(dt):
            form ='{:02}:{:02}:{:02}'.format( dt.hour, dt.minute, dt.second)
            return form
        
        sunset = self.timegrid[0]
        prepare = sunset - dt.timedelta(minutes=30)
        prepare_form = dt_to_form(prepare)
        sunset_form = dt_to_form(self.sunset)
        sunrise_form = dt_to_form(self.sunrise+ dt.timedelta(minutes = 20))
        flat_form = dt_to_form_flat(self.sunrise + dt.timedelta(minutes = 40))
        expset = list(set([exptime.split(',') for exptime in list(set(self.all_targetlist['exptime']))][0]))
        
        with open(filename, 'w') as f:
            
            # Preparation
            f.write('; Prepare for the observation\n')
            f.write('#WAITUNTIL 1, {}\n'.format(prepare_form))
            
            # Cooling
            f.write('\n; Cooling\n')
            f.write('#CHILL {}, {}\n'.format(ccdtemp, 1.0))
            
            # Calibration frames
            f.write('\n; Calibration frames\n')
            if dark:
                for exptime in expset:
                    f.write("#COUNT {}\n#INTERVAL {}\n#DARK\n".format(9, exptime)) # dark
            if bias:
                f.write("#COUNT {}\n#BIAS\n".format(9)) # bias
            
            # Image frames
            f.write('\n; Start of the evening\n')
            f.write('#WAITUNTIL 1, {}\n'.format(sunset_form))
            
            
            #f.write('#AFINTERVAL {}\n'.format(self.calib_interval//60)) # AF every 3 hrs
            
            f.write('\n; Targeting')

            total = 0
            timeslots = self.all_targetlist.group_by('group').groups
            for timeslot in timeslots:
                
                first_slottime = timeslot['exp_start(UT)'][0]
                first_slottime_form = dt_to_form(first_slottime)
                f.write('\n#WAITUNTIL 1, {}\n\n'.format(first_slottime_form))
                f.write('\n#AUTOFOCUS\n')
                for target in timeslot:
                    if target['priority'] == 5.0:
                        f.write('\n;------------------ToO Observation------------------\n')
                        name = target['name']
                        ra = target['ra']
                        dec = target['dec']
                        exptime = target['exptime']
                        filter_ = target['filter']
                        counts = target['counts']
                        binning = target['binning']
                        too_obstime = target['exp_start(UT)']
                        too_obstime_form = dt_to_form(too_obstime)
                        f.write('\n#WAITUNTIL 1, {}\n'.format(too_obstime_form))
                        #f.write('#POINTING\n')
                        f.write('#COUNT {}\n'.format(counts))
                        f.write('#INTERVAL {}\n'.format(exptime))
                        f.write('#BINNING {}\n'.format(binning))
                        f.write('#FILTER {}\n'.format(filter_))
                        f.write('{}\t{}\t{}\n'.format(name,ra,dec))
                        f.write('#NOPREVIEW \n')
                        f.write('#NOSOLVE\n')
                        f.write('#NOPOINTING\n')
                        f.write('\n;----------------ToO Observation end----------------\n')
                    else:
                        name = target['name']
                        ra = target['ra']
                        dec = target['dec']
                        exptime = target['exptime']
                        filter_ = target['filter']
                        counts = target['counts']
                        binning = target['binning']
                        #f.write('#POINTING\n')
                        f.write('#COUNT {}\n'.format(counts))
                        f.write('#INTERVAL {}\n'.format(exptime))
                        f.write('#BINNING {}\n'.format(binning))
                        f.write('#FILTER {}\n'.format(filter_))
                        f.write('{}\t{}\t{}\n'.format(name,ra,dec))
                        f.write('#NOPREVIEW \n')
                        f.write('#NOSOLVE\n')
                        #f.write('#NOPOINTING\n')
                        f.write('\n')
                    total += 1
            
            f.write('#QUITAT {}\n'.format(sunrise_form))
            if flat:
                f.write('#CHAIN autoflat_script.txt')
                with open(os.path.dirname(self.rtspath)+'/autoflat_script.txt', 'w') as f_flat:
                    f_flat.write('\n#WAITUNTIL 1, {}\n\n'.format(flat_form))
                    f_flat.write('#DAWNFLATS\n')
                    f_flat.write('\n; Closing\n')
                    f_flat.write('#SHUTDOWN\n')            
            else:
                f.write('\n; Closing\n')
                f.write('#SHUTDOWN\n')

        return filename

    def add_ToO(self,
                ID,
                ra,
                dec,
                obstime,
                exptime = '120,120',
                counts = '5,5',
                filter_ = 'g,r',
                binning = '1,1'):
        replaced_tbl = self.all_targetlist.copy()
        t_distance = np.abs(obstime - replaced_tbl['exp_start(UT)'])
        replace_idx = np.argmin(t_distance)

        replaced_tbl[replace_idx]['name'] = ID
        replaced_tbl[replace_idx]['ra'] = ra
        replaced_tbl[replace_idx]['dec'] = dec
        replaced_tbl[replace_idx]['exp_start(UT)'] = obstime
        replaced_tbl[replace_idx]['exptime'] = exptime
        replaced_tbl[replace_idx]['counts'] = counts
        replaced_tbl[replace_idx]['filter'] = filter_
        replaced_tbl[replace_idx]['binning'] = binning
        replaced_tbl[replace_idx]['priority'] = 5.0
        totexp = self.calc_totexp(replaced_tbl[replace_idx])
        remove_idx_candidate = replace_idx + 1
        while replaced_tbl[remove_idx_candidate]['exp_start(UT)'] < obstime + dt.timedelta(seconds = int(totexp)):
            replaced_tbl.remove_row(remove_idx_candidate)
            remove_idx_candidate += 1
        self.all_targetlist = replaced_tbl
        return replaced_tbl
    
    def cut_ToOscript(self, 
                      scriptpath):
        with open(scriptpath, 'r') as f:
            lines = f.readlines()
            idx = np.where(np.array(lines) == ';------------------ToO Observation------------------\n')[0]
            f.close()
        if len(idx) == 0:
            print('ToO Observation is not included in the script')
        else:
            idx = idx[0]
            new_lines = lines[idx:]
            with open(scriptpath, 'w') as newf:
                for new_line in new_lines:
                    newf.write(new_line)
                f.close()

# %%go
#%%
if __name__ == '__main__': 
    
    start = '2023/04/07'
    subsequent_days = 20
    
    start_datetime = datetime.datetime.strptime(start, '%Y/%m/%d')
    for day in range(subsequent_days):
        start_str = datetime.datetime.strftime(start_datetime, '%Y/%m/%d')
        rtspath = rtsmaker(start = start_str, end = None, save_path = '/home/hhchoi1022/ACPscript/', catpath = '../../Config/alltarget_KCT_220904.txt', altlimit = 30)
        scriptmaker = ACP_scriptmaker(rtspath, overhead_ratio = 1.1, select_mode = 'optimal', dark = True, bias = True, flat = False)
        start_datetime += datetime.timedelta(days = 1)
        
    
    
    
  # %%
