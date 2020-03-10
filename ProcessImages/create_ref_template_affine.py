# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 08:28:24 2015

@author: Kristof Govaerts

This script uses an affine registration algorithm to repeatedly register a 
set of images to an atlas or reference image, and uses a brain mask (in atlas space)
to find where the brain is in native space. Brain masks are initially 
dilated quite a lot to account for imperfect registration, but are iteratively 
reduced in size as the registration gets more accurate.
"""

import os
import nibabel as nib
import pylab as pl
import numpy as np
from datetime import datetime
import glob
import Tkinter, tkFileDialog
import subprocess
import scipy.ndimage as ndi
from skimage import filter

#Globals
N_ITERATIONS = 1

#functions
def RegAndResample(ref, atl, mode='A'):
    '''Registers [atl] into [ref] space and resamples [ref] into [atl] space.
    Inputs should be image paths to feed into FSL.
    Can do affine (mode='A') or rigid (mode='R') registration.'''
    omat = os.path.join(os.getcwd(), 'temp_omat.mat')    
    out = os.path.join(os.getcwd(), ref.split('.')[0] + mode + '.nii.gz') 
    if mode == 'A': #affine
        flirt = 'flirt -in %s -ref %s -out %s -omat %s -bins 128 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear -v' %(ref, atl, out, omat)  
    elif mode == 'R': #rigid
        flirt = 'flirt -in %s -ref %s -out %s -omat %s -bins 128 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear -v' %(ref, atl, out, omat)  
    
    subprocess.call(flirt, shell=True)
    #subprocess.call(resample, shell=True)
    os.remove(omat)
    return out
    
def RegAndResample_NRR(ref, atl):
    '''Registers [atl] into [ref] space and resamples [ref] into [atl] space.
    Inputs should be image paths to feed into FSL.
    Does nonrigid registration.'''
 
    out = os.path.join(os.getcwd(), ref.split('.')[0] + 'NRR' + '.nii.gz') 
    dfout = os.path.join(os.getcwd(), ref.split('.')[0] + 'DF' + '.nii.gz') 
    
    fnirt = "fnirt --ref=%s --in=%s --iout=%s --fout=%s --subsamp=4,2,2 --reffwhm=2,1,0, -infwhm=4,2,1" %(atl, ref, out, dfout)
    
    subprocess.call(fnirt, shell=True)
    #subprocess.call(resample, shell=True)
    return out
    
def otsu_filter(img, background=0):
    flat = (img.flatten())[img.flatten()>background]
    thresh = filter.threshold_otsu(flat)
    im2 = np.copy(img)
    im2[im2<=thresh]=0
    return im2
    
def reg_batch(refim, filelist, mode, zoom=1):
    '''Performs affine or rigid registration of 'filelist' to 'refim' and 
    resamples into 'refim' space. 'mode' should be 'rigid' or 'affine'.'''
    outlist = []
    for i, f in enumerate(filelist):
        print "file %s of %s." %(i+1, len(filelist))
        im = nib.load(f)
        data = im.get_data()
        affine = zoom_affine(im.get_affine(), zoom) 

        temp_nif = nib.Nifti1Image(im.get_data(), affine)
        nib.save(temp_nif, f[:-4] + 'RS.nii')
        out = os.path.join(rootdir, f[:-4]+'mode')  
        in_img = f[:-4] + 'RS.nii'
        if mode == 'NRR':
            out = RegAndResample_NRR(in_img, refim)
        else:
            out = RegAndResample(in_img, refim, mode=mode) #one registration step
        outlist.append(out)
        try:
            os.remove(in_img)
        except OSError:
            print "Could not remove %s. Remove manually." %in_img
            pass
    return outlist
    
def save_sum_img(filelist, filename):
    cim = nib.load(filelist[0]).get_data()
    sum_img = np.zeros((cim.shape[0], cim.shape[1], cim.shape[2], len(filelist)))
    
    for i, f in enumerate(filelist):
        cim = np.squeeze(nib.load(f).get_data())
        sum_img[:,:,:,i] = cim  
       
    nif1 = nib.Nifti1Image(sum_img, nib.load(filelist[0]).get_affine())
    nib.save(nif1, 'test.nii')       
     
    adata = (np.mean(sum_img, axis=-1)).astype('int32') 
    adata[adata<7000] = 0  
    nif = nib.Nifti1Image(adata, nib.load(filelist[0]).get_affine())
    nib.save(nif, filename)
    
def zoom_affine(affine, zoom):
    '''Multiplies the voxel size in an affine matrix by [zoom].'''
    aff=np.copy(affine)
    for i,l in enumerate(aff[:-1]):
        l[i]*=zoom
    return aff


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

print "Select your reference image (Will be used as float for the registration)"

root = Tkinter.Tk()
while True:
    refim = tkFileDialog.askopenfilename(initialdir=rootdir,title='Please select a file')
    root.destroy()
    break

#resave refim as .nii - better compatible with FSL
ri = nib.load(refim)
aff = ri.get_affine()

ni = nib.Nifti1Image(ri.get_data(), zoom_affine(ri.get_affine(), 10))
nib.save(ni, refim[:-4] + '_resave.nii')
refim = refim[:-4] + '_resave.nii'

#rootdir = '/media/sf_host/data/NRR/Raw_data/test_FSL'
#os.chdir(rootdir)
#refim = '/media/sf_host/data/NRR/Raw_data/Bl6_Ref_ER.hdr'
#maskim = '/media/sf_host/data/NRR/Raw_data/Bl6_ref_ER_BNL_lbl_BM.hdr'

filelist=glob.glob('*.hdr')
print "STEP 1: Performing rigid registration of all images to reference image."
rigid_list = reg_batch(refim, filelist, mode='R', zoom=10)
print "STEP 2: Averaging registered images to create rigid template."
rigid_temp = 'rigid_template.nii'
save_sum_img(rigid_list, rigid_temp)
print "STEP 3: Performing affine registration of all images to reference image."
affine_list = reg_batch(rigid_temp, filelist, mode='A', zoom=10)
print "STEP 4: Averaging registered images to create affine template."
affine_temp = 'affine_template.nii'
save_sum_img(affine_list, affine_temp)
print "STEP 5: Performing nonrigid registration of all images to affine reference image."
nrr_list = reg_batch(affine_temp, affine_list, mode='NRR', zoom=1)
print "STEP 6: Averaging registered images to create nonrigid template."
nrr_temp = 'nrr_template.nii'
save_sum_img(nrr_list, nrr_temp)

        
subprocess.call()