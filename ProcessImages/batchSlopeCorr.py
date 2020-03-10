# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 16:16:08 2014

@author: Kristof Govaerts
Slope correction and optional threshold for fMRI data.
"""

import scipy.ndimage as ndi
from skimage.filter import threshold_otsu
from imageTypes import *
import Tkinter, tkFileDialog
root = Tkinter.Tk()
import glob

THRESH = True #whether or not to remove background by thresholding the data. Should be False if supplying a mask
THRESH_fact = 0.8 #only adjust if you're missing parts of your img
MASK_EXT = False #if there is an associated brain mask for volume 1, False otherwise
MC = True #motion correction in FSL. ONLY WORKS IN LINUX BOX
ZOOM = 20 #False if no zoom necessary

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

for i, f in enumerate(filelist):
    print "File %s of %s: %s" %(i+1, len(filelist), f)
    ana = Bruker2AnalyzeImg(f[:-4], repetition_avg=False)
    ana.correctSlope() #use slope corr from imageTypes 
    
    affine = ana.affine
    #add a zoom factor

    if MASK_EXT != False:
        mim = np.squeeze(nib.load(f[:-4] + MASK_EXT).get_data())
        ana.pdata[:,:,:,0] *= mim
    
    if ZOOM != False:
        print "Adding zoom of factor %s." %ZOOM
        affine = ana.affine * np.diag([ZOOM,ZOOM,ZOOM,1.])
        affine[:,3] = affine[:,3]
        ana.affine=affine
    
    if THRESH == True:
        print "Applying threshold at %s * Otsu threshold per volume" %THRESH_fact
        for sl in range(ana.shape[-1]):
            mask = ana.pdata[:,:,:,sl]#for each rep seperately, registration comes after
            mask = ndi.gaussian_filter(mask, 1)
            t = threshold_otsu(mask)
            ana.pdata[:,:,:,sl][mask<THRESH_fact * t] = 0
    
    if MC == True:
        print "Motion correction in progress"
        nif = nib.Nifti1Image(ana.pdata, ana.affine)
        nib.save(nif, 'temp.nii')
        args = '-in %s' %os.path.join(os.getcwd(), 'temp.nii')   
        mc = 'mcflirt -in %s -out %s' %('temp.nii', f[:-4] + '_MC.nii.gz')   
        print mc
        subprocess.call(mc, shell=True)
        ana.pdata = nib.load(f[:-4] + '_MC.nii.gz').get_data()
        os.remove('temp.nii')
    
    if MASK_EXT != False:
        print "Applying supplied brain mask."
        for v in range(ana.shape[-1]):
            ana.pdata[:,:,:,v] *= mim
        
    nif = nib.Nifti1Image(ana.pdata, affine)
    
    nib.save(nif, f[:-4] + "_SC")
    del ana #clear memory, fMRI data is big