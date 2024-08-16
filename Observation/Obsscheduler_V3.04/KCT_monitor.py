#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 09:43:25 2022

@author: hhchoi10222
"""

import os, glob
import time
from astropy.io import fits
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
#%%
class Target:
    def __init__(self, 
                 watchdir = '/data6/obsdata/KCT_STX16803/zip', 
                 cpdir = '/data6/obsdata/KCT_STX16803',
                 check_period = 10,
                 recipients = ['hhchoi1022@gmail.com']):
        
        self.observer = Observer()
        self.watchdir = watchdir
        self.period = check_period
        self.cpdir = cpdir
        self.tmpdir = f'{watchdir}/temp'
        os.makedirs(self.tmpdir, exist_ok= True)
        self.recipients = recipients

    def run(self):
        event_handler = Handler(self.cpdir, self.recipients, self.period)
        self.observer.schedule(event_handler, self.watchdir, recursive=False)
        self.observer.start()
        try:
            while True:
                time.sleep(self.period)
        except:
            self.observer.stop()
            #self.observer.join()
            
class Handler(FileSystemEventHandler):
    
    def __init__(self, cpdir, recipients, check_period):    
        self.cpdir = cpdir
        self.recipients = recipients
        self.period = check_period

    def on_created(self, event):
        
        if (event.src_path).endswith('.zip'):
            # wait for fininshing copying file 
            historicalSize = -1
            while (historicalSize != os.path.getsize(event.src_path)):
                historicalSize = os.path.getsize(event.src_path)
                time.sleep(self.period)
            self.unzip(event.src_path, modify_header= True)
            
    def change_hdr(self, folderpath):
        filelist = glob.glob(f"{folderpath}/*")
        for filepath in filelist:            
            try:
                hdr = fits.getheader(filepath)
                if 'OBJECT' not in hdr.keys():
                    fits.setval(filepath, 'OBJECT', value = '')
                    print(f'{os.path.basename(filepath)} is modified')
            except:
                pass
    
    def unzip(self, event_path, modify_header : bool = True):
        import datetime
        os.system(f'chmod 777 {event_path}')
        os.system(f'unzip -d {self.tmpdir} {event_path}')
        foldername = os.path.basename(event_path).split('.')[0]
        folderpath = f'{self.tmpdir}/{foldername}'
        if modify_header:
            self.change_hdr(folderpath)
        os.system(f'mv {folderpath} {self.cpdir}')
        os.system(f'chmod 777 {folderpath}')
        # mailing
        subject = '[KCT] Observation completed'
        text = f'When : {datetime.datetime.today().strftime("%Y-%m-%d")}\nFrom : {event_path}\nTo : {folderpath}'
        self.send_mail(subject, text, self.recipients)
        
    def send_mail(self, subject, text, recipients):
    # Sending the mail
        import smtplib
        from email.mime.text import MIMEText
        from email.mime.multipart import MIMEMultipart
        
        smtp_server = smtplib.SMTP(host ='smtp.gmail.com', port = 587)
        smtp_server.ehlo()
        smtp_server.starttls()
        my_ID = 'hhchoi1022@snu.ac.kr'
        my_pw = 'evhivfptcuqkaant'
        smtp_server.login(user=my_ID,password=my_pw)
        
        body = text
        recipients_join = ','.join(recipients)
        my_msg = MIMEMultipart()
        txt1 = MIMEText('================ Observation completed ================\n')
        tbl = MIMEText(body)
        txt2 = MIMEText('\n================ Observation completed ================\n\n')
        txt3 = MIMEText('This message is automatically sent after the observation\n\n Hyeonho Choi')
        
        my_msg['Subject'] =subject
        my_msg['To'] = recipients_join
        my_msg.attach(txt1)
        my_msg.attach(tbl)
        my_msg.attach(txt2)
        my_msg.attach(txt3)

        smtp_server.sendmail(from_addr = my_ID,
                              to_addrs = recipients,
                              msg = my_msg.as_string())
        smtp_server.quit()
        
#%%
if __name__ == "__main__":
    watchdir = '/data6/obsdata/KCT_STX16803/zip'
    cpdir = '/data6/obsdata/KCT_STX16803'
    check_period = 3600
    recipients = ['hhchoi1022@gmail.com']
    w = Target(watchdir = watchdir, cpdir = cpdir, check_period = check_period, recipients = recipients)
    w.run()
    
