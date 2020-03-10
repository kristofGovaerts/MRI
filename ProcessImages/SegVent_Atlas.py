# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 08:28:24 2015

@author: Kristof Govaerts

This script iteratively registers all files in a directory you choose to a template image.
It creates a template image using rigid registration, then repeats the process with affine regisration,
updates the template and floating images, and repeats the process with nonrigid registration.
Adjust the N_ITERATIONS parameter to choose the number of images. NEUROMORPH_LOC needs to point
to the installation directory of Neuromorph (or any directory with all executables for NiftyReg).
THRESH can be adjusted if the background value for your images is not 0. 
"""

import os
import nibabel as nib
import numpy as np
import glob
import Tkinter, tkFileDialog
import subprocess
import shutil
import scipy.ndimage as ndi
from skimage.filter import threshold_otsu
from scipy.misc import imsave

#Globals
REFIM = 'C:/Users/u0091609/Desktop/BSI_normalized/normalized_all/BNLref_DS.hdr'
REFLBL = 'C:/Users/u0091609/Desktop/BSI_normalized/normalized_all/BNL_avgA_inv_DS.hdr'
VENT_IND = [10] #label index of ventricles - can be list of multiple labels
FILE_EXT = '*_mid_1*.nii'

NEUROMORPH_LOC = "C:/Program Files/NeuroMorph/Packages/MIC/Applications/Modules/Macros/KULNeuroMorph"
#default pars
BENDING_ENERGY = 0.05
GRID_SPACING = -5 #final grid spacing - if negative, is in voxels, if positive in mm
THRESH=0
FTYPE='float'

reg_ala = os.path.join(NEUROMORPH_LOC, 'reg_aladin.exe')
reg_f3d = os.path.join(NEUROMORPH_LOC, 'reg_f3d.exe')
reg_transform = os.path.join(NEUROMORPH_LOC, 'reg_transform.exe')
reg_resample = os.path.join(NEUROMORPH_LOC, 'reg_resample.exe')
reg_jacobian = os.path.join(NEUROMORPH_LOC, 'reg_jacobian.exe')
reg_average = os.path.join(NEUROMORPH_LOC, 'reg_average.exe')

#functions
def RegAndResample(flt, atl, filename, mode='A'):
    '''Registers [atl] into [flt] space and resamples [flt] into [atl] space.
    Inputs should be image paths.
    Can do affine (mode='A') or rigid (mode='R') registration.'''
    omat = filename[:-4]+'.txt'    
    out = filename 
    if mode == 'A': #affine    
        args = '-flo %s -ref %s -noSym' %(os.path.join(os.getcwd(), flt), atl)         
    elif mode == 'R': #rigid
        args = '-flo %s -ref %s -noSym -rigOnly' %(os.path.join(os.getcwd(), flt), atl)   
        
    flirt = [reg_ala] + args.split()    
 
    process = subprocess.call(flirt, shell=True)
    os.rename('outputResult.nii', out) #NiftyReg saves as outputResult by default
    os.rename('outputAffine.txt', omat)

    return out
    
def RegAndResample_NRR(flt, atl, ext='', ln=4, sx=-5, be=0.005, smooF=0):
    '''Registers [atl] into [flt] space and resamples [flt] into [atl] space.
    Inputs should be image paths.
    Does nonrigid registration.'''
 
    out = os.path.join(os.getcwd(), flt.split('.')[0] + 'NRR' + str(ext) +  '.nii') 
    dfout = os.path.join(os.getcwd(), flt.split('.')[0] + 'CPP' + str(ext) + '.nii') 
    
    args = '-flo %s -ref %s -ln %s -sx %s -be %s -smooF %s' %(os.path.join(os.getcwd(), flt), atl, str(ln), str(sx), str(be), str(smooF))     
    fnirt = [reg_f3d] + args.split()  
    
    subprocess.call(fnirt, shell=True)
    os.rename('outputResult.nii', out) #NiftyReg saves as outputResult by default
    os.rename('outputCPP.nii', dfout)
    return out, dfout
    
def otsu_filter(img, background=0):
    flat = (img.flatten())[img.flatten()>background]
    thresh = threshold_otsu(flat)
    im2 = np.copy(img)
    im2[im2<=thresh]=0
    return im2
    
def apply_transform(t, ref, flo, filename):
    '''Inputs:
    t: Filename of a control point image 
    ref: Reference image
    flo: Floating image to which to apply the transform
    filename: Name to give the transformed floating image.
    Output:
    Filename of the transformed floating image.'''
    args = '-ref %s -flo %s -trans %s -res %s -inter 0' %(ref, flo, t, filename)     
    tf = [reg_resample] + args.split()  
    subprocess.call(tf, shell=True)
    return filename    
    
def get_ventricles(im, lbl, vent_value):
    vmask = nib.load(lbl).get_data()
    vmask2 = np.zeros(vmask.shape)

    for v in vent_value:
        vmask2[vmask==v] = 1
    data = np.squeeze(nib.load(im).get_data())
    aff = nib.load(im).get_affine()
    
    vm_dil = ndi.morphology.binary_dilation(vmask2, iterations=12)
    vol_dil = np.sum(vm_dil)         
    
    of = otsu_filter(data * vm_dil) 
    of[of>0]=1    
    
    vol_otsu = np.sum(of)    
    if vol_otsu > vol_dil/3.0: #check if the otsu filter worked - if not, probably a whole block of ventricle was selected
        for i in range(11,3,-1):
            print "Could not find ventricles. Trying smaller dilation factor", i
            vm_dil = ndi.morphology.binary_dilation(vmask2, iterations=i)
            vol_dil = np.sum(vm_dil)    
            of = otsu_filter(data * vm_dil) 
            of[of>0]=1    
            vol_otsu = np.sum(of)   
            if vol_otsu < vol_dil/3.0:
                break
    if vol_otsu > vol_dil/3.0:
        print "Could not find ventricles in this image."
    
    nif = nib.Nifti1Image(of, aff)
    nib.save(nif, im[:-4] + '_Vent.nii')
    return im[:-4] + '_Vent.nii'
    
    
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

#rootdir = '/media/sf_host/data/NRR/Raw_data/test_FSL'
#os.chdir(rootdir)
#refim = '/media/sf_host/data/NRR/Raw_data/Bl6_Ref_ER.hdr'
#maskim = '/media/sf_host/data/NRR/Raw_data/Bl6_ref_ER_BNL_lbl_BM.hdr'

filelist=glob.glob(FILE_EXT)

for i, f in enumerate(filelist):
    print "File %s of %s." %(i+1, len(filelist))
    im = nib.load(f)
    nz = im.get_affine()
    nd = im.get_data()
    nd[nd<0]=0
    nif = nib.Nifti1Image(im.get_data(), nz)
    f2 = f.replace('.nii', 'Z.nii')
    nib.save(nif, f2)
    
    affim = RegAndResample(REFIM, f2, f[:-4]+'AffAtl.nii', mode='A')
    affmat = affim.replace('.nii', '.txt')
    rs_lbl = apply_transform(affmat, affim, REFLBL, f.replace('.nii', '_lbls.nii'))
    
    vents = get_ventricles(f2, rs_lbl, VENT_IND)
    vents = f2[:-4] + '_Vent.nii'
    vd = nib.load(vents).get_data()
    md = np.squeeze(nd) * vd
    png = np.mean(md, axis=-1) #MIP across Z-axis
    imsave(f[:-4] + '_mip.png', png)
    
    os.remove(affim)
    
    

