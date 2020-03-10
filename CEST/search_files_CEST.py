# -*- coding: utf-8 -*-
"""
Created on Tue May 29 17:10:51 2018

@author: Kristof

Searches all files
"""
import os
import pandas as pd
import Tkinter, tkFileDialog
root = Tkinter.Tk()

METH = '##$Method=<User:rare_cest>'

print "Enter the scan directory."

while True:
    rootdir = tkFileDialog.askdirectory(initialdir="/",title='Please select a directory')
    if os.path.isdir(rootdir) is True: #Checks if entered dir exists
        os.chdir(rootdir)
        root.destroy()
        break
    else:
        print "Pathname invalid. Try again."
        continue

os.chdir(rootdir)

def search_files(directory='.', extension=''):
    '''Gets all files as well as nested files.'''
    olist = []
    extension = extension.lower()
    for dirpath, dirnames, files in os.walk(directory):
        for name in files:
            if extension and name.lower().endswith(extension):
                olist.append(os.path.join(dirpath, name))
            elif not extension:
                olist.append(os.path.join(dirpath, name))
                
    return olist

def get_subdirectories():
    return next(os.walk('.'))[1]

def list_scans_with_meth(meth):
    c1= [os.path.join(os.getcwd(), i) for i in get_subdirectories()]
    ol =[]
    for d in c1:
        os.chdir(d)
        for f in get_subdirectories():
            try:
                with open(os.path.join(f, 'method')) as fi:
                    txt = fi.readlines()
            except IOError:
                continue
            if meth in ''.join(txt):
                ol.append(os.path.join(d, f))
        os.chdir('..')
    return ol

def output_cest_parfile():
    frame = pd.DataFrame()
    frame[0] = list_scans_with_meth(METH)
    frame[1] = [os.path.split(s)[0] + '_' + os.path.split(s)[1] for s in frame[0]]
    frame.columns=['path','out_name']
    frame.to_csv('CEST_parfile.csv', sep=';',index=False)
    
output_cest_parfile()