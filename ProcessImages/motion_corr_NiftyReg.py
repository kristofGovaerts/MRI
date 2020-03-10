# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 08:28:24 2015

@author: Kristof Govaerts

"""

import os
import nibabel as nib
import numpy as np
import glob
import Tkinter, tkFileDialog
import subprocess
import scipy.ndimage as ndi
from skimage.filters import gaussian_filter
from pandas import rolling_mean

#Globals
FUNC_EXT='7.nii'
THRESH=5.0
NEUROMORPH_LOC = "C:/Program Files/NeuroMorph/Packages/MIC/Applications/Modules/Macros/KULNeuroMorph"
SMOOTH = 0.75
T_SMOOTH = 5

reg_ala = os.path.join(NEUROMORPH_LOC, 'reg_aladin.exe')
reg_transform = os.path.join(NEUROMORPH_LOC, 'reg_transform.exe')
reg_resample = os.path.join(NEUROMORPH_LOC, 'reg_resample.exe')

#functions
def RegAndResample(flt, atl, mode='A'):
    '''Registers [atl] into [flt] space and resamples [flt] into [atl] space.
    Inputs should be image paths.
    Can do affine (mode='A') or rigid (mode='R') registration.'''
    omat = flt.split('.')[0] + mode + '.txt'    
    out = flt.split('.')[0] + mode + '.nii'
    if mode == 'A': #affine    
        args = '-flo %s -ref %s' %(flt, atl)         
    elif mode == 'R': #rigid
        args = '-flo %s -ref %s -rigOnly' %(flt, atl)   
        
    flirt = [reg_ala] + args.split()    
 
    process = subprocess.call(flirt, shell=True)
    os.rename('outputResult.nii', out) #NiftyReg saves as outputResult by default
    os.rename('outputAffine.txt', omat)

    return out, omat    

def temporal_smoothing(data, t_smooth):
    out_arr = np.array(data)
    for x in range(data.shape[0]):
        for y in range(data.shape[1]):
            for z in range(data.shape[2]):
                out_arr[x,y,z,:] = rolling_mean(out_arr[x,y,z,:], t_smooth, min_periods=1)   
    return out_arr

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
filelist=glob.glob('*' + FUNC_EXT)

print "STEP 1: Rigid reg of functional image to anatomical image."
for n, f in enumerate(filelist):
    print "File %s of %s." %(n+1, len(filelist))
    root = f[:-len(FUNC_EXT)+1]
    im = nib.load(f)
    aff = im.get_affine()
    data = im.get_data()
    if THRESH != 0:
        data[data<THRESH] = 0.0
    if SMOOTH != 0:
        for t in range(data.shape[-1]):
            for z in range(data.shape[-2]):
                data[:,:,z,t] = gaussian_filter(data[:,:,z,t], SMOOTH)
                
    #save ref img
    ref_img = np.mean(data, axis=-1)
    ref_name = root + '_av.nii'
    ref_nif = nib.Nifti1Image(ref_img, aff)
    nib.save(ref_nif, ref_name)
    
    #save source imgs:
    reglist = []
    matlist = []
    for t in range(data.shape[-1]):
        print "Volume %s of %s." %(t+1, data.shape[-1])
        num = "%03d" % (t)
        sourcename = root + num + '.nii'
        #save source
        sourcenif = nib.Nifti1Image(data[:,:,:,t], aff)
        nib.save(sourcenif, sourcename)

        #reg source
        out, omat = RegAndResample(sourcename, ref_name, mode='R')
        reglist += [out]
        matlist += [omat]
        os.remove(sourcename)
    out_arr = np.zeros(data.shape)
    for i, r in enumerate(reglist):
        out_arr[:,:,:,i] = nib.load(r).get_data()
        os.remove(r)
        os.remove(matlist[i])
        
    #temporal smoothing
    out_arr = temporal_smoothing(out_arr, T_SMOOTH)
     
     
    out_nif = nib.Nifti1Image(out_arr, aff)
    nib.save(out_nif, root + '_MC.nii')