#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 16:12:10 2022

@author: hhchoi1022
"""
import os
os.chdir('/home/hhchoi1022/Desktop/Gitrepo/Reference_make')
import imaplib
import imapclient
import email
from email.utils import parseaddr
from email.header import decode_header
import re
import datetime
import pandas as pd
from astropy.table import Table
from astropy.io import ascii
from Observation import cross_match
from Observation import to_skycoord
import astropy.units as u
from astropy.coordinates import SkyCoord


def AT_notifier(sepadist = 3600):

    # log-in
    imap_server = imaplib.IMAP4_SSL(host='imap.gmail.com',port='993')
    
    
    my_ID = 'hhchoi1022@snu.ac.kr'
    my_pw = 'evhivfptcuqkaant'
    imap_server.login(user=my_ID,password=my_pw)
    

    # check mail-box
    resp_list,mailbox_list = imap_server.list()

    # Set mailbox for parsing
    mailbox = "TNS"
    mailbox_code = imap_server.select(mailbox)
    
    code1, mail_all = imap_server.search(None,'All')
    
    # Parsing
    result = {}
    mail_ids = mail_all[0].split()
    
    def HtmltoText(data):
        data = re.sub(r'&nbsp;','',data)
        data = re.sub(r'</.*?>','\n',data)
        data = re.sub(r'<.*?>', '', data)
        return data
    
    mid = mail_ids[-1]
    code2, data = imap_server.fetch(mid,'(RFC822)')
    mail = {}
    msg = email.message_from_string(data[0][1].decode('utf-8'))

    mail['From'] = parseaddr(msg['From'])[1]
    mail['To'] = parseaddr(msg['To'])[1]
    mail['Date'] = datetime.datetime.strptime(msg['date'], '%a, %d %b %Y %H:%M:%S +0000')

    subject = decode_header(msg['Subject'])
    mail['Subject'] = subject[0][0]

    if msg.is_multipart():
        for part in msg.walk():
            if part.get_content_type() == 'text/plain':
                body = part.get_payload(decode=True)
                body = body.decode()    
    else:
        body = str(msg.get_payload(decode= True))
    
    mail['Body'] = str(HtmltoText(body))
    imap_server.close()

    # Constructing daily transients table
    

    tnslist = mail['Body'].split('Data source group')[:-1]
    ATtable = pd.DataFrame(index=range(0,len(tnslist)),columns=['TNSName','DetDate','DetMag','RA','Dec','reporter'])

    for i, tns in enumerate(tnslist):
        ra = re.findall('RA=(\d\d:\d\d:\d\d)', tns)[0]
        dec = re.findall('DEC=([+-]\d\d:\d\d:\d\d)',tns)[0]
        date = re.findall('Discovery date=(\d\d\d\d-\d\d-\d\d)',tns)[0]
        mag = re.findall('Discovery mag=(\d+)',tns)[0]
        name = re.findall('\d\d\d\d\w+',tns)[0]
        reporter = re.findall('Reporting group: (\w+)',tns)[0]
        
        ATtable['TNSName'][i] = name
        ATtable['DetDate'][i] = date
        ATtable['DetMag'][i] = mag
        ATtable['RA'][i] = ra
        ATtable['Dec'][i] = dec
        ATtable['reporter'][i] = reporter

    

    # Matching with alltargets
    ATtable = Table.from_pandas(ATtable)
    
    alltargets = ascii.read('/data1/obsroutine/alltarget.dat')
    alltargets = alltargets[(alltargets['priority']!=1.0)|(alltargets['priority']!=1.1)|(alltargets['priority']!=1.5)]
    targetlist = alltargets['obj']
    SKY_coord = SkyCoord(alltargets['ra'], alltargets['dec'],unit=(u.hourangle,u.deg))
    AT_coord = SkyCoord(ATtable['RA'], ATtable['Dec'],unit = (u.hourangle,u.deg))
    matched_skyidx, matched_atidx, no_match_idx = cross_match(SKY_coord,AT_coord, sepadist)
    matchedtarget = alltargets[matched_skyidx]
    matchedat = ATtable[matched_atidx]
    radectarget = SkyCoord(matchedtarget['ra'],matchedtarget['dec'],unit=(u.hourangle,u.deg))
    radecat = SkyCoord(matchedat['RA'],matchedat['Dec'],unit=(u.hourangle,u.deg))
    sepdist = SkyCoord.separation(radectarget,radecat).arcsecond.round(0)
    matchedat['host'] = matchedtarget['obj']
    matchedat['separation[arcsec]'] = sepdist
    return matchedat


#%%

def send_mail(recipients,htmlfile='/data1/tmp/tns.html' ):
    # Sending the mail
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.base import MIMEBase
    from email.mime.multipart import MIMEMultipart
    from email import encoders
    import os
    import codecs
    
    smtp_server = smtplib.SMTP(host ='smtp.gmail.com', port = 587)
    smtp_server.ehlo()
    smtp_server.starttls()
    my_ID = 'hhchoi1022@snu.ac.kr'
    my_pw = 'evhivfptcuqkaant'
    smtp_server.login(user=my_ID,password=my_pw)
    
    subject = '[IMSNG FIELD] TNS Transient Alert'
    f = codecs.open(htmlfile,'r')
    body = f.read()
    
    
    recipients_join = ','.join(recipients)
    
    my_msg = MIMEMultipart()
    tbl = MIMEText(body,'html')
    txt = MIMEText('This is automatically made AT list from TNS daily digest','plain')
    my_msg['Subject'] =subject
    my_msg['To'] = recipients_join
    my_msg.attach(tbl)
    my_msg.attach(txt)
    
    smtp_server.sendmail(from_addr = my_ID,
                          to_addrs = recipients,
                          msg = my_msg.as_string())
    smtp_server.quit()
    print(f'Mail sent to {recipients_join}')

#%%
recipients = ['hhchoi1022@gmail.com']
table = AT_notifier(3600)
table.write('/data1/tmp/tns.html',format = 'ascii.html',overwrite=True)
table
if len(table) !=0:
    send_mail(recipients)




# Sending mail
# Automation
