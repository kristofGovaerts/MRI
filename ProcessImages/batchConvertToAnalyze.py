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

SMOOTHING = 1.5 #Gaussian smoothing sigma value = FWHM. 
MODE = "2D" #2D if smoothing on a slice-by-slice basis, 3D if smoothing per volume

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
        
    nif = nib.AnalyzeImage(data, affine)
    
    nib.save(nif, f[:-4] )    

    del img #clear memory