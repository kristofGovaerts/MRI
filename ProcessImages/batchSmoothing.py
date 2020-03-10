# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 16:16:08 2014

@author: u0091609

Batch smoothing tool for fMRI data.
"""

import pylab as pl
import scipy.ndimage as ndi
import nibabel as nib
import Tkinter, tkFileDialog
root = Tkinter.Tk()
import glob
import os
import numpy as np

SMOOTHING = [1, 1, 2.5] #Gaussian smoothing sigma value = FWHM. 
MODE = "3D" #2D if smoothing on a slice-by-slice basis, 3D if smoothing per volume

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

filelist = glob.glob("*.nii")

for i, f in enumerate(filelist):
    print "File %s of %s: %s" %(i+1, len(filelist), f)
    img = nib.load(f)
    affine = img.get_affine()
    data = img.get_data()
    if len(data.shape) == 3:
        data = np.reshape(data, list(data.shape) + [1]) #reshape so that you have a 4D array 
    mask = np.array(data)
    mask[mask>0] = 1
        
    if MODE == "2D":
        print "Performing 2D smoothing with a kernel of %s" %SMOOTHING
        for i in range(data.shape[-1]):
            for j in range(data.shape[-2]):                
                data[:,:,j,i] = ndi.filters.gaussian_filter(data[:,:,j,i], SMOOTHING)
                
    elif MODE == "3D":
        print "Performing 3D smoothing with a kernel of %s" %SMOOTHING
        for i in range(data.shape[-1]):              
            data[:,:,:,i] = ndi.filters.gaussian_filter(data[:,:,:,i], SMOOTHING)
            
    data *= mask
        
    nif = nib.Nifti1Image(data, affine)
    
    nib.save(nif, f[:-4] + '_S' + str(SMOOTHING).replace('.','p'))    

    del img #clear memory, fMRI data is big
