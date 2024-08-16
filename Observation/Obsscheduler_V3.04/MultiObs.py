#%%
from obsscheduler import ObsScheduler
from astropy.io import ascii
import json
from astropy.time import Time
from googlesheet import GoogleSheet
from NINA_Observation import NINA_observation
from tqdm import tqdm
import re
from astropy.table import Table, vstack
#%%
def multi_obs_spec(name_sheet = '231102_multi',
                  name_project : str = 'GECKO',
                  specmode_configfolder = './config/7DT/specmode/',
                  specmode_filename = 'SpecMode_6.config',
                  ccdtemp : float = -10,
                  gain : int = 2750,
                  delay_minute : float = 10,
                  safe : bool = True,
                  sunalt : bool = False,
                  wait_until : bool = False,
                  bias_before_obs : bool = True,
                  bias_after_obs : bool = True,
                  bias_count : int = 9,
                  dark_before_obs : bool = False,
                  dark_after_obs : bool = False,
                  dark_count : int = 9,        
                  autofocus : bool = False,
                  autofocus_interval_in_hour : bool = False,
                  autofocus_at_init : bool = True,
                  autofocus_filterchange : bool = True,
                  autofocus_every_minute : bool = False,
                  autofocus_every_minute_value : int = 60,
                  entire_night : bool = False,
                  obsdate : Time = None,
                  savepath = './NINAscript/',
                  update_ggsheet : bool = True,
                  show_altaz_each_target : bool = False
                  ):
    """
    name_sheet = 'S240422ed_UPDATE_ToO'
    name_project : str = 'GECKO'
    specmode_configfolder = './config/7DT/specmode/'
    specmode_filename = 'SpecMode_10.config'
    ccdtemp : float = -15
    delay_minute : float = 10
    safe : bool = True
    sunalt : bool = False
    wait_until : bool = False
    bias : bool = True
    dark_before_obs : bool = True
    dark_after_obs : bool = True
    autofocus : bool = False
    autofocus_interval_in_hour : bool = False
    autofocus_at_init : bool = False
    autofocus_filterchange : bool = True
    autofocus_every_minute : bool = False
    autofocus_every_minute_value : int = 60
    entire_night : bool = True
    obsdate : Time = None
    savepath = './NINAscript/'
    update_ggsheet : bool = True
    """
    sheet = GoogleSheet()
    tbl_multi = sheet.get_sheet_data(name_sheet, format_ = 'Table')
    
    if obsdate == None:
        obsdate_strlist = re.findall("(\d{2})(\d{2})(\d{2})", name_sheet)[0]
        obsdate_str = re.findall("(\d{2}\d{2}\d{2})", name_sheet)[0]
        obsdate = Time(f'20{obsdate_strlist[0]}-{obsdate_strlist[1]}-{obsdate_strlist[2]}')
    else:
        obsdate_str = '%.4d%.2d%.2d'%(obsdate.datetime.year, obsdate.datetime.month, obsdate.datetime.day)
    
    specmode_filepath = specmode_configfolder + specmode_filename
    with open(specmode_filepath, 'r') as f:
        specmode_pathdict = json.load(f)
    specmode_keys = list(specmode_pathdict.keys())
    specmode_paths = list(specmode_pathdict.values())
    
    # Define the specmode 
    specmode_dict = dict()
    for mode in specmode_keys:
        specmode_configfile = specmode_configfolder + specmode_pathdict[mode]
        with open(specmode_configfile, 'r') as f:
            specmode = json.load(f)
        specmode_dict[mode] = specmode
    name_telescopes = list(specmode_dict[next(iter(specmode_dict))].keys())
    
    # Target table for each telescope
    target_tbl_dict = dict()
    for name_telescope in name_telescopes:
        target_tbl = tbl_multi.copy()
        filter_all = []
        for i, target in enumerate(target_tbl):
            mode = target['mode']
            if mode not in specmode_keys:
                filters_str = mode
            else:
                filters = specmode_dict[mode][name_telescope]
                filters_str = ','.join(filters)
            filter_all.append(filters_str)
        target_tbl['filter'] = filter_all
        target_tbl_dict[name_telescope] = target_tbl
    
    # ObsScheduler for each telescope
    schedule_dict = dict()
    for name_telescope in tqdm(name_telescopes):
        target_tbl = target_tbl_dict[name_telescope]
        scheduler = ObsScheduler(target_db = target_tbl, date = obsdate, name_project = name_project, name_telescope = name_telescope, entire_night = entire_night)
        schedule_dict[name_telescope] = scheduler.scheduler(delay_minute = delay_minute, autofocus = autofocus, autofocus_at_init = autofocus_at_init, autofocus_interval_in_hour = autofocus_interval_in_hour)
        
    for name_telescope in tqdm(name_telescopes):
        all_tbl = schedule_dict[name_telescope].all
        schedule_tbl = schedule_dict[name_telescope].scheduled
        if len(schedule_tbl) == 0:
            raise ValueError("No scheduled targets found")
        target_tbl = schedule_tbl[schedule_tbl['obj']!='autofocus']
        obs_exptime = [exptime.split(',') for exptime in target_tbl['exptime']]
        exptime_all = list(set([value for sublist in obs_exptime for value in sublist]))
        obs_filter = [filter_.split(',') for filter_ in target_tbl['filter']]
        filter_all = list(set([value for sublist in obs_filter for value in sublist]))
        obs_binning = [binning.split(',') for binning in target_tbl['binning']]
        binning_all = list(set([value for sublist in obs_binning for value in sublist]))


        # Make script based on the scheduled_dict =================================================
        # =========================================================================================
            
        nina = NINA_observation(name_telescope = name_telescope, gain = gain)
        nina.initialize_container(safe = safe, sunalt = sunalt)
        nina.wait_sun_alt(18, name_container = 'Wait')
        nina.cool(ccdtemp, duration_minute = 5, name_container = 'Wait')
        nina.unpark()
        
        # Preparation
        nina.wait_sun_alt(-10, name_container = 'Preparation')
        if bias_before_obs:
            for binning in binning_all:
                nina.bias(count = bias_count, binning = binning)
        if dark_before_obs:
            for binning in binning_all:
                for exptime in exptime_all:
                    nina.dark(count = dark_count, exptime = exptime, binning = binning)
                    pass
        nina.wait_sun_alt(-15, name_container = 'Preparation')
        
        # Observation
        if str(schedule_tbl['filter'][0]).split(',')[0] == 'g':
            nina.slewRADec(ra = float(schedule_tbl['ra'][0]), dec = float(schedule_tbl['dec'][0]), name_container= 'Sequential')
            nina.switch_flt('g', name_container = 'Sequential')
            nina.autofocus(name_container = 'Sequential')
        nina.target(target_tbl = schedule_tbl, wait_until = wait_until, autofocus_filterchange = autofocus_filterchange, autofocus_every_minute = autofocus_every_minute, autofocus_every_minute_value = autofocus_every_minute_value)
        
        if bias_after_obs:
            for binning in binning_all:
                nina.bias(count = bias_count, binning = binning, name_container= 'Sequential')
        if dark_after_obs:
            for binning in binning_all:
                for exptime in exptime_all:
                    nina.dark(count = dark_count, exptime = exptime, binning = binning, name_container = 'Sequential')
        
        # Shutdown (No cooler off)
        nina.slewAltAz(alt = 40, az = 300, name_container = 'Sequential')
        nina.tracking_off(name_container = 'Sequential')
        nina.write(savepath = savepath, filename = f'{obsdate_str}-{name_telescope}-Spec.json')
        # Make script based on the scheduled_dict =================================================
        # =========================================================================================

    if show_altaz_each_target:
        from mainobserver import mainObserver
        from maintarget import mainTarget
        for target in schedule_dict['7DT_01'].scheduled:
            targ = mainTarget(name_telescope= '7DT_01', name_project= name_project, observer = mainObserver(name_telescope = '7DT_01', name_project = name_project), target_ra = float(target['ra']), target_dec = float(target['dec']), target_name = target['obj'])
            targ.staralt(utctime = target['exec_start'])

    if update_ggsheet:
        for column in tbl_multi.colnames:
            tbl_multi[column] = all_tbl[column]
        for column in ['moonsep', 'scheduled', 'maxalt', 'fraction_obs', 'transit']:
            tbl_multi[column] = all_tbl[column]
        for column in ['transit']:
            tbl_multi[column] = [str(value) for value in all_tbl[column]]
        for column in ['filter']:
            tbl_multi[column] = ''
        #tbl_multi.sort('ra')
        tbl_multi['ra'] = tbl_multi['ra'].astype('float')
        tbl_multi['weight'] = tbl_multi['weight'].astype('float')
        tbl_multi.sort(['weight', 'ra'])
        
        sheet.write_sheet_data(sheet_name = name_sheet, data = tbl_multi, append = False)
            
        

#%%
import astropy.units as u
A = multi_obs_spec(name_sheet = '240526_multi',
                    name_project  = 'GECKO',
                    specmode_configfolder = './config/7DT/specmode/',
                    specmode_filename = 'SpecMode_10.config',
                    ccdtemp  = -10,
                    gain = 2750,
                    delay_minute  = 10,
                    wait_until  = False,
                    safe = False,
                    sunalt = False,
                    bias_before_obs = False,
                    bias_after_obs = True,
                    bias_count = 18,
                    dark_before_obs = False,
                    dark_after_obs = True,
                    dark_count = 18,                           
                    autofocus  = False,
                    autofocus_at_init  = False,
                    autofocus_filterchange  = True,
                    autofocus_every_minute  = False,
                    autofocus_every_minute_value  = 30,
                    entire_night  = True,
                    obsdate = Time.now(),
                    savepath = './NINAscript/',
                    update_ggsheet = False,
                    show_altaz_each_target = True)


#%%
def multi_obs_search(name_sheet = '231129_multi',
                     name_telescopes = ['7DT_01','7DT_02','7DT_03','7DT_05','7DT_06','7DT_07','7DT_08','7DT_09','7DT_10,','7DT_11'],
                     name_project : str = 'Tile',
                     ccdtemp : float = -10,
                     gain = 2750,
                     delay_minute : float = 0,
                     safe : bool = True,
                     sunalt : bool = False,
                     wait_until : bool = False,
                     bias_before_obs : bool = False,
                     bias_after_obs: bool = True,
                     bias_count : bool = 18,
                     dark_before_obs : bool = False,
                     dark_after_obs : bool = True,
                     dark_count : bool = 18,
                     autofocus : bool = False,
                     autofocus_interval_in_hour : bool = False,
                     autofocus_at_init : bool = False,
                     autofocus_filterchange : bool = True,
                     autofocus_every_minute : bool = True,
                     autofocus_every_minute_value : int = 45,
                     entire_night : bool = False,
                     obsdate : Time = None,
                     savepath = './NINAscript/',
                     update_ggsheet : bool = True,
                     show_altaz_each_target : bool = False
                     ):
    """
    name_sheet = 'S240422ed_ToO'
    name_telescopes = ['7DT_01','7DT_02','7DT_03','7DT_05','7DT_06','7DT_07','7DT_08','7DT_09','7DT_10','7DT_11']
    name_project : str = 'Tile'
    ccdtemp : float = -10
    gain = 2750
    delay_minute : float = 0
    safe : bool = True
    sunalt : bool = False
    wait_until : bool = False
    bias_before_obs : bool = False
    bias_after_obs: bool = True
    bias_count : bool = 18
    dark_before_obs : bool = False
    dark_after_obs : bool = True
    dark_count : bool = 18
    autofocus : bool = False
    autofocus_interval_in_hour : bool = False
    autofocus_at_init : bool = False
    autofocus_filterchange : bool = True
    autofocus_every_minute : bool = True
    autofocus_every_minute_value : int = 45
    entire_night : bool = False
    obsdate : Time = None
    savepath = './NINAscript/'
    update_ggsheet : bool = True
    show_altaz_each_target : bool = False
    """
    sheet = GoogleSheet()
    tbl_search = sheet.get_sheet_data(name_sheet, format_ = 'Table')
    
    if obsdate == None:
        obsdate_strlist = re.findall("(\d{2})(\d{2})(\d{2})", name_sheet)[0]
        obsdate_str = re.findall("(\d{2}\d{2}\d{2})", name_sheet)[0]
        obsdate = Time(f'20{obsdate_strlist[0]}-{obsdate_strlist[1]}-{obsdate_strlist[2]}')
    else:
        obsdate_str = '%.4d%.2d%.2d'%(obsdate.datetime.year, obsdate.datetime.month, obsdate.datetime.day)
    
    target_tbls = [tbl_search[i::len(name_telescopes)] for i in range(len(name_telescopes))]
    target_tbl_dict = dict()
    for i, name_telescope in enumerate(name_telescopes):
        target_tbl = target_tbls[i]
        target_tbl['filter'] = target_tbl['mode']
        target_tbl['unit_scope'] = name_telescope
        target_tbl_dict[name_telescope] = target_tbl

    # ObsScheduler for each telescope
    schedule_dict = dict()
    for name_telescope in tqdm(name_telescopes):
        target_tbl = target_tbl_dict[name_telescope]
        if len(target_tbl) > 0:
            scheduler = ObsScheduler(target_db = target_tbl, date = obsdate, name_project = name_project, name_telescope = name_telescope, entire_night = entire_night)
            schedule_dict[name_telescope] = scheduler.scheduler(delay_minute = delay_minute, autofocus = autofocus, autofocus_at_init = autofocus_at_init, autofocus_interval_in_hour = autofocus_interval_in_hour)
    
    update_tbl = Table()
    for name_telescope in tqdm(name_telescopes):
        all_tbl = schedule_dict[name_telescope].all
        update_tbl = vstack([all_tbl, update_tbl])
        schedule_tbl = schedule_dict[name_telescope].scheduled
        if len(schedule_tbl) == 0:
            raise ValueError("No scheduled targets found")
        target_tbl = schedule_tbl[schedule_tbl['obj']!='autofocus']
        obs_exptime = [exptime.split(',') for exptime in target_tbl['exptime']]
        exptime_all = list(set([value for sublist in obs_exptime for value in sublist]))
        obs_filter = [filter_.split(',') for filter_ in target_tbl['filter']]
        filter_all = list(set([value for sublist in obs_filter for value in sublist]))
        obs_binning = [binning.split(',') for binning in target_tbl['binning']]
        binning_all = list(set([value for sublist in obs_binning for value in sublist]))


        # Make script based on the scheduled_dict =================================================
        # =========================================================================================
            
        nina = NINA_observation(name_telescope = name_telescope, gain = gain)
        nina.initialize_container(safe = safe, sunalt = sunalt)
        nina.wait_sun_alt(0, name_container = 'Wait')
        nina.cool(ccdtemp, duration_minute = 5, name_container = 'Wait')
        nina.unpark()

        # Preparation
        nina.wait_sun_alt(-10, name_container = 'Preparation')
        if bias_before_obs:
            for binning in binning_all:
                nina.bias(count = bias_count, binning = binning)
        if dark_before_obs:
            for binning in binning_all:
                for exptime in exptime_all:
                    nina.dark(count = dark_count, exptime = exptime, binning = binning)
                    pass
        nina.target(target_tbl = schedule_tbl, wait_until = wait_until, autofocus_filterchange = autofocus_filterchange, autofocus_every_minute = autofocus_every_minute, autofocus_every_minute_value = autofocus_every_minute_value)
        
        if bias_after_obs:
            for binning in binning_all:
                nina.bias(count = bias_count, binning = binning, name_container= 'Sequential')
        if dark_after_obs:
            for binning in binning_all:
                for exptime in exptime_all:
                    nina.dark(count = dark_count, exptime = exptime, binning = binning, name_container = 'Sequential')
        
        #nina.slewAltAz(alt = 40, az = 300) 
        nina.write(savepath = savepath, filename = f'{obsdate_str}-{name_telescope}-Spec.json')
    
        # Make script based on the scheduled_dict =================================================
        # =========================================================================================

    if show_altaz_each_target:
        from mainobserver import mainObserver
        from maintarget import mainTarget
        for name_telescope in name_telescopes:
            for target in schedule_dict[name_telescope].scheduled:
                targ = mainTarget(name_telescope= name_telescope, name_project= name_project, observer = mainObserver(name_telescope = '7DT_01', name_project = name_project), target_ra = float(target['ra']), target_dec = float(target['dec']), target_name = target['obj'])
                targ.staralt(utctime = target['exec_start'])

    if update_ggsheet:
        update_tbl['weight'] = [float(value) for value in update_tbl['weight']]
        update_tbl.sort('weight')
        for column in ['unit_scope']:
            tbl_search[column] = update_tbl[column]
        for column in tbl_search.colnames:
            tbl_search[column] = update_tbl[column]
        for column in ['moonsep', 'scheduled', 'maxalt', 'fraction_obs', 'transit']:
            tbl_search[column] = update_tbl[column]
        for column in ['transit']:
            tbl_search[column] = [str(value) for value in update_tbl[column]]
        for column in ['filter']:
            tbl_search[column] = ''
        tbl_search['weight'] = tbl_search['weight'].astype(float)
        tbl_search.sort(['weight', 'ra'])
        sheet.write_sheet_data(sheet_name = name_sheet, data = tbl_search, append = False)
            

#%%
A = multi_obs_search(name_sheet = 'S240422ed_ToO',
                     name_telescopes = ['7DT_01','7DT_02','7DT_03','7DT_05','7DT_06','7DT_07','7DT_08','7DT_09','7DT_10','7DT_11'],
                     name_project  = 'GECKO',
                     ccdtemp  = -10,
                     gain = 2750,
                     delay_minute  = 0,
                     safe  = True,
                     sunalt  = False,
                     wait_until  = False,
                     bias_before_obs  = False,
                     bias_after_obs = True,
                     bias_count  = 18,
                     dark_before_obs  = False,
                     dark_after_obs  = True,
                     dark_count  = 18,
                     autofocus  = False,
                     autofocus_interval_in_hour  = False,
                     autofocus_at_init  = False,
                     autofocus_filterchange  = True,
                     autofocus_every_minute  = True,
                     autofocus_every_minute_value  = 45,
                     entire_night  = False,
                     obsdate  = None,
                     savepath = './NINAscript/',
                     update_ggsheet  = True,
                     show_altaz_each_target  = False)
# %%
