#%%
from obsscheduler import ObsScheduler, ScriptMaker
from googlesheet import GoogleSheet
from astropy.io import ascii
import json
from astropy.time import Time
from googlesheet import GoogleSheet
from NINA_Observation import NINA_observation
from tqdm import tqdm
from tqdm import trange
import re
from astropy.table import Table, vstack

class focus_test(NINA_observation):
    def __init__(self,
                 base_name = 'Base',
                 filtpath = './config/7DT/filtinfo',
                 name_telescope : str = '7DT_01'):
        super().__init__(base_name = base_name, filtpath = filtpath, name_telescope = name_telescope)
        self.name_telescope = name_telescope
        
    def sequence(self,
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
                autofocus_at_init : bool = True,
                focus_offset_value = 500,
                focus_offset_step = 4,
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
        else:
            for info_key in ['exptime', 'filter_', 'count', 'binning']:
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
        self.move_focus_relative(focus_offset = -focus_offset_value * focus_offset_step//2, name_container= target_name)
        for i in range(focus_offset_step+1):
            
            for i in range(len(target_info['filter_'])):
                ft = target_info['filter_'][i]
                et = target_info['exptime'][i]
                ct = target_info['count'][i]
                bt = target_info['binning'][i]
                sel = self.scriptmaker.se_maker(self.filtpath, ft, ct, ct, et, 1, 0, bt, bt, 'LIGHT')
                self.script = self.scriptmaker.add_action(sel, target_name, self.script)        
            self.move_focus_relative(focus_offset= focus_offset_value, name_container= target_name)
    
    def run(self,
            target_name : str = 'FocusTest',
            ra_h = 6,
            ra_m = 46,
            ra_s = 0,
            dec_d = -5,
            dec_m = 27,
            dec_s = 0,
            exptime = 10,
            filter_ = 'r',
            count = 1,
            binning = 1,
            **kwargs
            ):
        self.initialize_container()
        self.cool(ccdtemp, duration_minute = 5, name_container = 'Wait')
        self.unpark()
        self.sequence(target_name = target_name, ra_h = ra_h, ra_m= ra_m, ra_s=ra_s, dec_d= dec_d, dec_m= dec_m, dec_s= dec_s, exptime = exptime, filter_ = filter_, count = count, binning= binning, **kwargs)
        self.write(savepath = './NINAscript/', filename = f'Focustest.json')
#%%
f = focus_test()
# %%
ccdtemp = -15
bias = False
dark_before_obs = False
ra_h = 7
ra_m = 24
ra_s = 30
dec_d = -31
dec_m = 51
dec_s = 0
exptime = 10
filter_ = 'r'
binning = 1
count = 1
f.run(target_name = 'FocusTest', ra_h = ra_h, ra_m= ra_m, ra_s=ra_s, dec_d= dec_d, dec_m= dec_m, dec_s= dec_s, exptime = exptime, filter_ = filter_, count = count, binning= binning, focus_offset_value = 500, focus_offset_step = 8)
# %%
