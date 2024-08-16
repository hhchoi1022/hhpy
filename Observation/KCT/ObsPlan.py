# ACP Observing plan generator

import os
import sys
import glob
from astropy.io import ascii
from astropy.time import Time

class GenPlan :
    def KCT(input_name='/data1/target/KCT/2021/IMSNG-20210311-rts_vis-KCT.txt'):
        """
        time = 1 : Starting script including #CHILL #AUTOFOCUS #UNPARK(?)
        time = 2 : No special function
        time = 3 : Finishing script including #CHILLOFF #PARK or #TRACKOFF
        
        ;#AUTOFOCUS
        #POINTING
        #NOPREVIEW
        #COUNT 3,3,3
        #INTERVAL 120,120,120
        #BINNING 1,1,1
        #FILTER r,g,i 
        NGC1566	04:20:0.398	-54:56:16.120
        """
        rtscat  = ascii.read(input_name)
        project = input_name.split('/')[-1].split('-')[0]
        date    = input_name.split('/')[-1].split('-')[1]
        num = 1

        # Observing routines
        band    = 'B,V,R'
        exp     = '60,60,60'
        count   = '3,3,3'
        binning = '1,1,1'

        # target
        name = rtscat['name']
        ra   = rtscat['ra']
        dec  = rtscat['dec']

        # special target
        special1 ='M90'

        f = open('{}-{}{}.txt'.format(date, project, num), 'w+')
        f.write('#AUTOFOCUS\n')
        f.write('\n')
        for n,r,d in zip(name, ra, dec): 
            f.write('#POINTING\n')
            f.write('#COUNT {}\n'.format(count))
            f.write('#INTERVAL {}\n'.format(exp))
            f.write('#BINNING {}\n'.format(binning))
            f.write('#FILTER {}\n'.format(band))
            f.write('{}\t{}\t{}\n'.format(n,r,d))
            f.write('\n')
        f.close()
    def SAO(rtscat_name, date='', project='', num='', catalog=True):
        from astropy.time import TimezoneInfo  # Specifies a timezone
        from datetime import datetime   
        import pytz

        if catalog == False:
            rtscat  = ascii.read(rtscat_name)
            project = rtscat_name.split('/')[-1].split('-')[0]
            date    = rtscat_name.split('/')[-1].split('-')[1]
        elif catalog == True :
            rtscat = rtscat_name
        #num = 1
        # Observing routines
        band    = 'R'
        exp     = '30'
        count   = '12'
        binning = '1'
    
        # target
        name = rtscat['name']
        ra   = rtscat['ra']
        dec  = rtscat['dec']
        rise = rtscat['rise(LT)']
        
        utc = TimezoneInfo()    # Defaults to UTC
        utc_plus_9h = TimezoneInfo(utc_offset=9*u.hour)  # UTC+9
        local = pytz.timezone("Asia/Seoul")

        f = open('{}-{}-SAO_KL4040_HDR{}.txt'.format(date, project, num), 'w+')
        f.write('#AUTOFOCUS\n')
        f.write('\n')
        for n,r,d,rt in zip(name, ra, dec, rise): 

            naive = datetime.strptime(date[0:4]+'-'+date[4:6]+'-'+date[6:8]+' '+rt+':00', "%Y-%m-%d %H:%M:%S")
            local_dt = local.localize(naive, is_dst=None)
            utc_dt = local_dt.astimezone(pytz.utc)
            utc_dt_fmt = utc_dt.strftime("%Y-%m-%d %H:%M:%S")

            f.write('#WAITUNTIL 1, {}\n'.format(utc_dt_fmt.split(' ')[-1]))
            f.write('#NOPOINTING\n')
            f.write('#COUNT {}\n'.format(count))
            f.write('#INTERVAL {}\n'.format(exp))
            f.write('#BINNING {}\n'.format(binning))
            f.write('#FILTER {}\n'.format(band))
            f.write('{}\t{}\t{}\n'.format(n,r,d))
            f.write('\n')
        f.close()    

    def SAO_year(rtslist_name='/data1/target/SAO/IMSNG/2022/IMSNG*.txt'):
        """
        1. 20h < ra < 00h, 00h < ra < 07h
        2. 07h < ra < 12.5h
        3. 12.5h < ra < 20h
        """
        from astropy.time import Time  
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        rtslist = glob.glob(rtslist_name); rtslist.sort()
        for rts in rtslist:
            
            project = rts.split('/')[-1].split('-')[0]
            date    = rts.split('/')[-1].split('-')[1]
            incat = ascii.read(rts)
            ra, dec = incat['ra'], incat['dec']
            incat['radeg']  = incat['priority']
            incat['decdeg'] = incat['priority']
            coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            for i in range(len(coord)) :
                incat['radeg'][i] = (coord[i]).ra.value
                incat['decdeg'][i] = (coord[i]).dec.value
    
            p1 = incat[ ((incat['radeg'] > 20*15) & (incat['radeg'] < 24*15)) | ((incat['radeg'] > 0*15) & (incat['radeg'] < 7*15))]
            p2 = incat[(incat['radeg'] > 7*15) & (incat['radeg'] < 12.5*15)]
            p3 = incat[(incat['radeg'] > 12.5*15) & (incat['radeg'] < 20*15)]
            GenPlan.SAO(p1, date=date, project='IMSNG', num=1, catalog=True)
            GenPlan.SAO(p2, date=date, project='IMSNG',num=2, catalog=True)
            GenPlan.SAO(p3, date=date, project='IMSNG',num=3, catalog=True)
            #GenPlan.SAO(rts, catalog=False)
    
    
    def GenPlanAll(input_name='/home/lim9/anaconda3/lib/python3.7/site-packages/lgpy/data/alltarget_KCT.dat'):
        """
        ;#AUTOFOCUS9
        #POINTING
        #NOPREVIEW
        #COUNT 5,5,5
        #INTERVAL 60,60,60
        #BINNING 1,1,1
        #FILTER B,V,R
        NGC1566	04:20:0.398	-54:56:16.120
    
        # Scheduler script examples
        ; ACPX Horizon = 40
        ; ACPX moonavoid=120,14
        ; ACPS HoUrAnGlE = -3.5 , 4.0
        ;   ACPX Horizon  = none
        ; ACPS Priority = 10
        """
        rtscat  = ascii.read(input_name)
        #project = input_name.split('/')[-1].split('-')[0]
        #date    = input_name.split('/')[-1].split('-')[1]
        num = 1
    
        # Observing routines
        """
        band         = 'B,V,R'
        exp          = '60,60,60'
        count        = '5,5,5'
        binning      = '1,1,1'
        Preview      = False
        MoonAvoid    = '40,140' # Distance, Width --> 140 makes all moon distance == 40 deg
        """
        band         = 'g,r'
        exp          = '120,120'
        count        = '5,5'
        binning      = '1,1'
        Preview      =  False
        MoonAvoid    = '40,14' # Distance, Width --> 140 makes all moon distance == 40 deg    
    
        # target
        name = rtscat['obj']
        ra   = rtscat['ra']
        dec  = rtscat['dec']
        pr   = rtscat['priority']
    
        pr[pr == 1.0] = 3
        pr[pr == 1.1] = 2
        pr[pr == 1.5] = 1
    
        f = open('IMSNG-KCT_STX16803-pr1_1.5-20211216.txt', 'w+')
        #f.write('#AUTOFOCUS\n')
        #f.write('\n')
        for n,r,d,p in zip(name, ra, dec, pr): 
            f.write('; ACPX moonavoid={}\n'.format(MoonAvoid))
            f.write('; ACPX Priority = {}\n'.format(round(p)))
            if Preview == False :
                f.write('#NOPREVIEW\n')
            f.write('#POINTING\n')
            f.write('#COUNT {}\n'.format(count))
            f.write('#INTERVAL {}\n'.format(exp))
            f.write('#BINNING {}\n'.format(binning))
            f.write('#FILTER {}\n'.format(band))
            f.write('{}\t{}\t{}\n'.format(n,r,d))
            f.write('\n')
        f.close()
        
    def GenPlanAll():
        band         = 'g,r'
        exp          = '120,120'
        count        = '10,10'
        binning      = '1,1'
        Preview      =  False
        #MoonAvoid    = '40,14' # Distance, Width --> 140 makes all moon distance == 40 deg    
    
        # target
        name = 'GRB210919A'
        ra   = '05:20:55'
        dec  = '01:16:27'
        pr   = 4
        f = open('IMSNG-KCT_STX16803-GRB210919A.txt', 'w+')
        
        #f.write('; ACPX moonavoid={}\n'.format(MoonAvoid))
        f.write('; ACPX Priority = {}\n'.format(round(pr)))
        if Preview == False :
            f.write('#NOPREVIEW\n')
        f.write('#POINTING\n')
        f.write('#COUNT {}\n'.format(count))
        f.write('#INTERVAL {}\n'.format(exp))
        f.write('#BINNING {}\n'.format(binning))
        f.write('#FILTER {}\n'.format(band))
        f.write('{}\t{}\t{}\n'.format(name,ra,dec))
        f.write('\n')
        f.close()
    
