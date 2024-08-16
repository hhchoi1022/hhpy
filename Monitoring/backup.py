#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 19:49:04 2022

@author: hhchoi1022
"""


from datetime import datetime
import os

def gettoday():
    today = datetime.strftime(datetime.today(), '%y%m%d')
    return today

def execute(command):
    os.system(command)
    print('Execution complete')

def backup(targetdir = '/home/hhchoi1022/Desktop/Gitrepo/', backupdir = '/data1/Backup/'):
    today = gettoday()
    backupdir_f = backupdir+today 
    os.makedirs(backupdir_f, exist_ok = True)
    command_cp = f'cp -r {targetdir}/* {backupdir_f}'
    execute(command_cp)
    print(f'Backup complete Target directory : {targetdir}\nBackup directory : {backupdir}{today}\n')


backup()
