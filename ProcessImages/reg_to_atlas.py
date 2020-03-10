# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 08:28:24 2015

@author: Kristof Govaerts

Batch script to register a group of images to an atlas. 
"""

import os
import nibabel as nib
import numpy as np
import glob
import Tkinter, tkFileDialog
import subprocess

#Globals
IM_EXT='_av.nii'
M_EXT = '_M.nii'
ATLAS = 'TV_Gy0_NC_pN1_4_1_SC_av.nii'
MODE = 'A'
AFFDIRECT = False
NOSYM = True

NEUROMORPH_LOC = "C:/Program Files/NeuroMorph/Packages/MIC/Applications/Modules/Macros/KULNeuroMorph"
#default pars
THRESH=7000
FTYPE='float'

reg_ala = os.path.join(NEUROMORPH_LOC, 'reg_aladin.exe')
reg_transform = os.path.join(NEUROMORPH_LOC, 'reg_transform.exe')
reg_resample = os.path.join(NEUROMORPH_LOC, 'reg_resample.exe')

#functions
def RegAndResample(flt, atl, mask=None, mode='A'):
    '''Registers [atl] into [flt] space and resamples [flt] into [atl] space.
    Inputs should be image paths.
    Can do affine (mode='A') or rigid (mode='R') registration.'''
    omat = flt.split('.')[0] + mode + '.txt'    
    out = flt.split('.')[0] + mode + '.nii'
    if mode == 'A': #affine    
        args = '-flo %s -ref %s' %(flt, atl)         
    elif mode == 'R': #rigid
        args = '-flo %s -ref %s -rigOnly' %(flt, atl)   
    
    if mask != None:
        args += ' -fmask ' + mask
    if AFFDIRECT:
        args += ' -affDirect'
    if NOSYM: 
        args += ' -noSym'
        
    flirt = [reg_ala] + args.split()    
    print flirt
 
    process = subprocess.call(flirt, shell=True)
    os.rename('outputResult.nii', out) #NiftyReg saves as outputResult by default
    os.rename('outputAffine.txt', omat)

    return out, omat

def apply_transform(t, ref, flo, filename, inter=1):
    '''Inputs:
    t: Filename of a control point image or affine transformation matrix
    ref: Reference image
    flo: Floating image to which to apply the transform
    filename: Name to give the transformed floating image.
    Output:
    Filename of the transformed floating image.'''
    args = '-ref %s -flo %s -trans %s -res %s -inter %s' %(ref, flo, t, filename, inter)     
    tf = [reg_resample] + args.split()  
    print tf
    subprocess.call(tf, shell=True)
    return filename        
    
def save_sum_img(filelist, filename, thresh=7000, ftype='int32'):
    cim = nib.load(filelist[0]).get_data()
    sum_img = np.zeros((cim.shape[0], cim.shape[1], cim.shape[2], len(filelist)))
    
    for i, f in enumerate(filelist):
        cim = np.squeeze(nib.load(f).get_data())
        sum_img[:,:,:,i] = cim  
       
    nif1 = nib.Nifti1Image(sum_img, nib.load(filelist[0]).get_affine())
    nib.save(nif1, 'test.nii')       
     
    adata = (np.mean(sum_img, axis=-1)).astype(ftype) 
    adata[adata<thresh] = 0 #for intensity normalized images, the background is usually not 0. 
    nif = nib.Nifti1Image(adata, nib.load(filelist[0]).get_affine())
    nib.save(nif, filename)
    
def resave_with_affine(fn, aff, fn2):
    '''Resaves image with a new affine matrix.'''
    im = nib.load(fn)
    nif = nib.Nifti1Image(im.get_data(), aff)
    nib.save(nif, fn2)
    
    
#Main loop
print "Enter the scan directory. This directory should contain all scans except for the reference image, and no other scans."
root = Tkinter.Tk()
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

filelist=glob.glob('*' + IM_EXT)
masklist = [f.replace(IM_EXT, M_EXT) for f in filelist]

outlist = []
omatlist = []

for i, f in enumerate(filelist):
    print "File %s of %s." %(i+1, len(filelist))
    mask = masklist[i]
    out, omat = RegAndResample(f, ATLAS, mask=mask, mode=MODE)
    bmask = apply_transform(omat, ATLAS, mask, f.replace(IM_EXT, '_BMR.nii'), inter=0)
    outlist += out
    omatlist += omat