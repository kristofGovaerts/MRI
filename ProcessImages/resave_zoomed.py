# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 08:28:24 2015

@author: Kristof Govaerts

This script is for fMRI registration. You need three files per scan:
1. The actual, multi-volume fMRI scan
2. The expanded EPI image, using the same settings as for the fMRI but with 
higher SNR (more averages) and a more slices (whole-brain). Should be bias field corrected
3. A high-resolution, non-distorted anatomical image

And also one reference image to register all images to.

 NEUROMORPH_LOC needs to point
to the installation directory of Neuromorph (or any directory with all executables for NiftyReg).
THRESH can be adjusted if the background value for your images is not 0. 
"""

import os
import nibabel as nib
import numpy as np
import glob
import Tkinter, tkFileDialog
root = Tkinter.Tk()

ZOOM_FACTOR = [20, 20, 20, 1]
    
def resave_with_affine(fn, aff, fn2):
    '''Resaves image with a new affine matrix.'''
    im = nib.load(fn)
    nif = nib.Nifti1Image(im.get_data(), aff)
    nib.save(nif, fn2)
    
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
filelist = glob.glob("*.hdr")
if len(filelist) == 0:
    filelist = glob.glob("*.nii")

for n, f in enumerate(filelist):
    print "File %s of %s." %(n+1, len(filelist))
    affine = nib.load(f).get_affine()
    n_affine = affine * np.diag(ZOOM_FACTOR)
    resave_with_affine(f, n_affine, f[:-4] + "Z.nii")
    