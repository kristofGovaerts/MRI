# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 16:16:08 2014

@author: u0091609

Batch averaging tool for fMRI data. For registration purposes.
"""
import nibabel as nib
import Tkinter, tkFileDialog
root = Tkinter.Tk()
import glob
import numpy as np
import os

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

filelist = glob.glob("*SC.nii")

for i, f in enumerate(filelist):
    print "File %s of %s: %s" %(i+1, len(filelist), f)
    img = nib.load(f)
    affine = img.get_affine()
    data = img.get_data()

    data= np.average(data, axis=-1)

    nif = nib.Nifti1Image(data, affine)
    
    nib.save(nif, f[:-4] + '_avg')    

    del img #clear memory, fMRI data is big
