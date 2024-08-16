    
#%%
from NINAconverter import NINA
import json
from astropy.coordinates import SkyCoord
from astropy.time import Time
from mainobserver import mainObserver
from NINAconverter import NINA
from googlesheet import GoogleSheet
from astropy.io import ascii
from obsscheduler import ObsScheduler
import astropy.units as u

class NINA_observation(NINA):
    def __init__(self,
                 base_name = 'Base',
                 filtpath = './config/7DT/filtinfo',
                 name_telescope : str = '7DT_01', 
                 gain : int = 2750
                 ):
        super().__init__()
        self.basename = base_name
        self.filtpath = filtpath
        self.nametelescope = name_telescope
        self.action = None
        self.trigger = None
        self.scriptmaker = None
        self.condition = None
        self.container = None
        self.filters = None
        self.script = None
        self.observer = None
        self.wait_first = True
        self.gain = gain
        
        self._set_actions()
        self._set_filtfile()
        self._set_observer()
        self.initialize_container()
    
    def _set_observer(self):
        self.observer = mainObserver(name_telescope = self.nametelescope, name_project = 'GECKO')
        
    def _set_actions(self):
        nina = NINA()
        self.action = nina.Action()
        self.trigger = nina.Trigger()
        self.scriptmaker = nina.Scriptmaker()
        self.condition = nina.Condition()
        self.container = nina.Container()
        
    def _set_filtfile(self):
        self.filtpath = self.filtpath + f'/{self.nametelescope}'
        filtfile = self.filtpath + '/filter.json'
        with open(filtfile, 'r') as f:
            json_data = json.load(f)
        self.filters = [json_data['Items']['$values'][1]['Items']['$values'][i]['Items']['$values'][0]['Filter']['_name'] for i in range(len(json_data['Items']['$values'][1]['Items']['$values']))]
            
    def initialize_container(self,
                             meridianflip : bool = True,
                             sunalt : bool = True,
                             sunalt_value : float = -18,
                             safe : bool = True,
                             ):
        cont_base = self.container.base(name = 'Base')
        cont_wait = self.container.sequentialContainer('Wait')
        cont_prep = self.container.sequentialContainer('Preparation')
        cont_seq = self.container.sequentialContainer('Sequential')
        self.script = self.scriptmaker.add_container(cont_wait,'Targets',cont_base)
        self.script = self.scriptmaker.add_container(cont_prep,'Targets',cont_base)
        self.script = self.scriptmaker.add_container(cont_seq,'Targets',cont_base)
        
        ####### Add condition here!
        if meridianflip:
            meridianflip = self.trigger.meridianflip()
            self.script = self.scriptmaker.add_trigger(meridianflip, 'Sequential', self.script)
        if sunalt:
            sunalt = self.condition.sun_alt(0, 0, 0, False, 0, 0, 0, sunalt_value, 3)
            self.script = self.scriptmaker.add_condition(sunalt, 'Sequential', self.script)
        if safe:
            safe = self.condition.safe()
            self.script = self.scriptmaker.add_condition(safe, 'Sequential', self.script)

    def wait_sun_alt(self, altitude = 0, name_container : str = 'Wait'):
        wtsun = self.action.wait_sunaltitude(0, 0, 0, False, 0, 0, 0, altitude, 3)
        self.script = self.scriptmaker.add_action(wtsun, name_container, self.script)
        
    def wait(self, ut, name_container : str = 'Sequential'):
        def ut_2_lt(ut):
            return self.observer.localtime(ut)
        lt = ut_2_lt(ut)
        if self.wait_first:
            wt = self.action.wait_time(lt.hour, lt.minute, 0, lt.second)
            self.wait_first = False
        else:
            wt = self.action.wait_time_dup(lt.hour, lt.minute, 0, lt.second)
        self.script = self.scriptmaker.add_action(wt, name_container, self.script)
    
    def tracking_on(self, name_container : str = 'Sequential'):
        track = self.action.set_track_side()
        self.script = self.scriptmaker.add_action(track, name_container, self.script)
        
    def tracking_off(self, name_container : str = 'Sequential'):
        track = self.action.set_track_stop()
        self.script = self.scriptmaker.add_action(track, name_container, self.script)  
    
    def cool(self, temperature : int = -10, duration_minute : int = 1, name_container : str = 'Preparation'):
        cool = self.action.camera_cooling(temperature, duration_minute)
        self.script = self.scriptmaker.add_action(cool, name_container, self.script)
    
    def warm(self, duration_minute : int = 3, name_container : str = 'End'):
        warm = self.action.camera_warming(duration_minute)
        self.script = self.scriptmaker.add_action(warm, name_container, self.script)

    def switch_flt(self, filt : str, name_container : str = 'Sequential'):
        switch = self.action.switch_filt(filt = filt)
        self.script = self.scriptmaker.add_action(switch, name_container, self.script)
    
    def park(self, name_container : str = 'End'):
        park = self.action.park()
        self.script = self.scriptmaker.add_action(park, name_container, self.script)
    
    def unpark(self, name_container : str = 'Preparation'):
        unpark = self.action.unpark()
        self.script = self.scriptmaker.add_action(unpark, name_container, self.script)
    
    def autofocus(self, name_container = 'Sequential'):
        autofocus = self.action.autofocus()
        self.script = self.scriptmaker.add_action(autofocus, name_container, self.script)
    
    def move_focus_relative(self, focus_offset : int, name_container = 'Sequential'):
        focus = self.action.move_focus_relative(focus_offset)
        self.script = self.scriptmaker.add_action(focus, name_container, self.script)
    
    def slewAltAz(self, alt : float, az : float, name_container = 'Sequential'):
        alt_deg = int(alt)
        alt_minute = int((alt - alt_deg)*60)
        alt_second = int((alt - alt_deg - alt_minute / 60) * 3600)
        az_deg = int(az)
        az_minute = int((az - az_deg)*60)
        az_second = int((az - az_deg - az_minute/60) * 3600)
        slew_aa = self.action.slew_AA(AzD_t = az_deg, AzM_t = az_minute, AzS_t = az_second, AltD_t = alt_deg, AltM_t = alt_minute, AltS_t = alt_second)
        self.script = self.scriptmaker.add_action(slew_aa, name_container, self.script)
        
    def slewRADec(self, ra : float, dec : float, name_container = 'Sequential'):
        skycoord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
        ra_hour = skycoord.ra.hms.h
        ra_minute = skycoord.ra.hms.m
        ra_second = skycoord.ra.hms.s
        
        dec_deg = skycoord.dec.dms.d
        dec_minute = skycoord.dec.dms.m
        dec_second = skycoord.dec.dms.s
        dec_negative = False
        slew_rd = self.action.slew_RD(RAH_t= ra_hour, RAM_t= ra_minute, RAS_t= ra_second, NegativeDec_t= dec_negative, DecD_t= dec_deg, DecM_t= dec_minute, DecS_t= dec_second)
        self.script = self.scriptmaker.add_action(slew_rd, name_container, self.script)
        
    def bias(self, count : int = 5, exptime : float = 0, binning : int = 1, name_container : str = 'Preparation'):
        seb = self.scriptmaker.se_maker(self.filtpath ,self.filters[0], count, count, exptime,self.gain,0,binning, binning,'BIAS')
        self.script = self.scriptmaker.add_action(seb, name_container, self.script)
    
    def dark(self, count : int = 5, exptime : float = 120.0, binning : int = 1, name_container : str = 'Preparation'):
        sed = self.scriptmaker.se_maker(self.filtpath ,self.filters[0], count, count, exptime,self.gain,0,binning, binning,'DARK')
        self.script = self.scriptmaker.add_action(sed, name_container, self.script)
        
    def flat(self, count : int = 5, exptime : float = 120.0, filter_ = 'g', binning : int = 1, azimuth = 300, altitude = 40, name_container : str = 'Sequential'):
        azi_deg = int(azimuth)
        azi_minute = int((azimuth - azi_deg)*60)
        azi_second = int(((azimuth - azi_deg)*60 - azi_minute) * 60)
        alt_deg = int(altitude)
        alt_minute = int((altitude - alt_deg)*60)
        alt_second = int(((altitude - alt_deg)*60 - alt_minute) * 60)
        slew = self.action.slew_AA(AzD_t= azi_deg, AzM_t = azi_minute, AzS_t = azi_second, AltD_t = alt_deg, AltM_t = alt_minute, AltS_t = alt_second)
        self.script = self.scriptmaker.add_action(slew, name_container, self.script)
        
        if not filter_ in self.filters:
            raise ValueError("")
        sef = self.scriptmaker.se_maker(self.filtpath , filter_, count, count, exptime,self.gain,0,binning, binning,'FLAT')
        self.script = self.scriptmaker.add_action(sef, name_container, self.script)
    
    def add_target(self,
                   target_name : str,
                   ra_h,
                   ra_m,
                   ra_s,
                   dec_d,
                   dec_m,
                   dec_s,
                   exptime,
                   filter_,
                   count,
                   binning,
                   name_container : str = 'Sequential',
                   autofocus_filterchange : bool = False,
                   autofocus_every_minute : bool = False,
                   autofocus_every_minute_value : int = 60,
                   autofocus_at_init : bool = False,
                   ):
        target_info = dict()
        target_info['exptime'] = exptime
        target_info['filter_'] = filter_
        target_info['count'] = count
        target_info['binning'] = binning
        if ',' in target_info['filter_']:
            target_info['filter_'] = target_info['filter_'].split(',')
            for info_key in ['exptime', 'count', 'binning']:
                try:
                    target_info[info_key] = target_info[info_key].split(',')
                    if len(target_info[info_key]) != len(target_info['filter_']):
                        target_info[info_key] = [target_info[info_key] for i in range(len(target_info['filter_']))]    
                except:
                    target_info[info_key] = [target_info[info_key] for i in range(len(target_info['filter_']))]

        if int(dec_d) < 0:
            dec_negative = True
        else:
            dec_negative = False

        DSO = self.container.DSO(target_name,0,ra_h,ra_m,ra_s,dec_d,dec_m,dec_s,dec_negative)
        self.script = self.scriptmaker.add_container(DSO, name_container, self.script)
        slew = self.action.slew_RD(RAH_t= ra_h, RAM_t= ra_m , RAS_t= ra_s, NegativeDec_t = dec_negative, DecD_t= dec_d, DecM_t= dec_m, DecS_t= dec_s)
        self.script = self.scriptmaker.add_action(slew, target_name, self.script)
        if autofocus_at_init:
            self.switch_flt(filt = target_info['filter_'][0], name_container= target_name)
            self.autofocus(name_container= target_name)
        if autofocus_filterchange:
            aff = self.trigger.AF_filterchan()
            self.script = self.scriptmaker.add_trigger(aff, target_name, self.script)
        if autofocus_every_minute:
            aft = self.trigger.AF_time(amt = autofocus_every_minute_value)
            self.script = self.scriptmaker.add_trigger(aft, target_name, self.script)

        for i in range(len(target_info['filter_'])):
            ft = target_info['filter_'][i]
            et = target_info['exptime'][i]
            ct = target_info['count'][i]
            bt = target_info['binning'][i]
            sel = self.scriptmaker.se_maker(self.filtpath, ft, ct, ct, et, self.gain, 0, bt, bt, 'LIGHT')
            self.script = self.scriptmaker.add_action(sel, target_name, self.script)        
    
    def target(self, target_tbl,
               autofocus_filterchange : bool = True,
               autofocus_every_minute : bool = False,
               autofocus_every_minute_value : int = 60,
               wait_until : bool = False):
        dup_idx = 1

        def ut_2_lt(ut):
            return self.observer.localtime(ut)

        for i, target in enumerate(target_tbl):

            previous_target_tbl = target_tbl[:i]
            if wait_until:
                ut_obs_start = Time(target['exec_start']).datetime
                self.wait(ut_obs_start, name_container= 'Sequential')
            if target['obj'] == 'autofocus':
                if i+1 != len(target_tbl):
                    target_next = target_tbl[i+1]
                else:
                    target_next = target_tbl[i-1]
                name = target_next['obj']
                ra = target_next['ra_hms'].split(':')
                dec = target_next['dec_dms'].split(':')
                exptime = target_next['exptime'].split(',')
                filter_ = target_next['filter'].split(',')
                count = target_next['count'].split(',')
                binning = target_next['binning'].split(',')
                if int(dec[0]) < 0:
                    dec_negative = True
                else:
                    dec_negative = False

                slew = self.action.slew_RD(RAH_t= ra[0], RAM_t= ra[1], RAS_t= ra[2], NegativeDec_t= dec_negative, DecD_t= dec[0], DecM_t= dec[1], DecS_t= dec[2])
                self.script = self.scriptmaker.add_action(slew, 'Sequential', self.script)
                self.autofocus()
            else:
  
                name = target['obj']
                if name in previous_target_tbl['obj']:
                    name += f'dup_{dup_idx}'
                    dup_idx += 1
                ra = target['ra_hms'].split(':')
                dec = target['dec_dms'].split(':')
                exptime = target['exptime'].split(',')
                filter_ = target['filter'].split(',')
                count = target['count'].split(',')
                binning = target['binning'].split(',')
                if int(dec[0]) < 0:
                    dec_negative = True
                else:
                    dec_negative = False
                DSO = self.container.DSO(name,0,ra[0],ra[1],ra[2],dec[0],dec[1],dec[2],dec_negative)
                self.script = self.scriptmaker.add_container(DSO, 'Sequential', self.script)
                slew = self.action.slew_RD(RAH_t= ra[0], RAM_t= ra[1] , RAS_t= ra[2], NegativeDec_t = dec_negative, DecD_t= dec[0], DecM_t= dec[1],DecS_t= dec[2])
                self.script = self.scriptmaker.add_action(slew, name, self.script)
                if autofocus_filterchange:
                    aff = self.trigger.AF_filterchan()
                    self.script = self.scriptmaker.add_trigger(aff, name, self.script)
                if autofocus_every_minute:
                    aft = self.trigger.AF_time(amt = autofocus_every_minute_value)
                    self.script = self.scriptmaker.add_trigger(aft, name, self.script)
                for i in range(len(filter_)):                        
                    ft = filter_[i]
                    et = exptime[i]
                    ct = count[i]
                    bt = binning[i]
                    sel = self.scriptmaker.se_maker(self.filtpath, ft, ct, ct, et, self.gain, 0, bt, bt, 'LIGHT')
                    self.script = self.scriptmaker.add_action(sel, name, self.script)
    
    def write(self,
              filename : str,
              savepath = './NINAscript/'
              ):
        self.scriptmaker.update_dict_ids(self.script)
        self.scriptmaker.provider_update(self.script)
        self.scriptmaker.remove_duplicate_targetname(self.script)
        filename = filename#f'/NINA_{self.nametelescope}.json'
        self.scriptmaker.write(savepath, filename, self.script)
        print(f"{filename} is saved")
    
    
    
    
    
        
    

#%%
if __name__ == '__init__':
    sheet = GoogleSheet()
    date = Time.now() + 0 * u.day
    target_tbl = sheet.get_sheet_data('format_IMSNG', format_ = 'Table')
    nina = NINA_observation()
    obsscheduler = ObsScheduler(target_db = target_tbl, date = date, name_project='GECKO', entire_night = True)
    schedule = obsscheduler.scheduler(autofocus_at_init = True, autofocus = True, duplicate_when_empty = True)

    # Wait
    nina.wait_sun_alt(0, name_container = 'Wait')

    # Preparation
    nina.wait_sun_alt(-10, name_container = 'Preparation')
    nina.cool()
    nina.unpark()
    nina.bias(count = 9)
    nina.dark(count = 9, exptime = 120)
    nina.wait((Time(schedule.scheduled[0]['exec_start'])).datetime, name_container = 'Preparation')

    # Sequential
    nina.target(schedule.scheduled, wait_until = True)
    nina.warm()
    nina.park()
    nina.write(f'IMSNG_{date.datetime.year}{date.datetime.month}{date.datetime.day}.json')

# %%
