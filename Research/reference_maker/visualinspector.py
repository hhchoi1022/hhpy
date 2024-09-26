#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:05:42 2021

@author: hhchoi10222
"""

import tkinter
from tkinter import filedialog
import glob, os
from multiprocessing import Process

target = 'ESO182-G010'
filterlist = ['g','r','i']

def remove_file(path):
    file_path = 'open'
    while type(file_path)== str:
        root = tkinter.Tk()
        root.withdraw()
        file_path = filedialog.askopenfilenames(parent=root,initialdir=path)
        for file in file_path:
            #os.system(f'rm{file}')
            print(f'{os.path.basename(file)} removed!')
            
for filter_ in filterlist:
    imkey = f'/data3/hhchoi1022/IMSNG/KCT_STX16803/selected/{target}/{filter_}/*.fits'
    imagelist = glob.glob(imkey)
    path = os.path.dirname(imagelist[0])
    imlist_sliced = [imagelist[i:i+40] for i in range(0,len(imagelist),40)]
    t = Process(target = remove_file,args = (path,))
    t.start()
    for imlist in imlist_sliced:
        command = 'ds9'
        for imname in imlist:
            command = command+' '+imname
        os.system(command)
        stop = input('stop?')
        if stop == '1':
            break