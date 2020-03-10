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
VENT_IND = [1] #label index of ventricles - can be list of multiple labels
FILE_EXT = '*_BrainD_M4_FSN1.hdr'
LBL_EXT = '_brain_M4_FS_TregS_bl6_cba_ref5it_NR_VL2.nii'

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
def otsu_filter(img, background=0):
    flat = (img.flatten())[img.flatten()>background]
    thresh = threshold_otsu(flat)
    im2 = np.copy(img)
    im2[im2<=thresh]=0
    return im2
    
def get_ventricles(im, lbl, vent_value):
    vmask = np.squeeze(nib.load(lbl).get_data())
    
    vmask2 = np.zeros(vmask.shape)
    
    for v in vent_value:
        print "Found label at index", v
        vmask2[vmask==v] = 1
    data = np.squeeze(nib.load(im).get_data())
    aff = nib.load(im).get_affine()
    
    vm_dil = ndi.morphology.binary_dilation(vmask2, iterations=4)
    vol_dil = np.sum(vm_dil)      
    
    of = otsu_filter(data * vm_dil) 
    of[of>0]=1    
    
    vol_otsu = np.sum(of) 
    if vol_otsu > vol_dil/1.2: #check if the otsu filter worked - if not, probably a whole block of ventricle was selected
        for i in range(3,2,-1):
            print "Could not find ventricles. Trying smaller dilation factor", i
            vm_dil = ndi.morphology.binary_dilation(vmask2, iterations=i)
            vol_dil = np.sum(vm_dil)    
            of = otsu_filter(data * vm_dil) 
            of[of>0]=1    
            vol_otsu = np.sum(of)   
            if vol_otsu < vol_dil/1.2:
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
    lbl = f.replace(FILE_EXT[1:], LBL_EXT)
    
    try:
        vents = get_ventricles(f, lbl, VENT_IND)
    except IOError:
        print "No label image found. Continuing."
        continue
    vents = f[:-4] + '_Vent.nii'
    vd = nib.load(vents).get_data()
    md = np.squeeze(nib.load(f).get_data()) * vd
    png = np.mean(md, axis=-1) #MIP across Z-axis
    imsave(f[:-4] + '_mip.png', png)
    
    

