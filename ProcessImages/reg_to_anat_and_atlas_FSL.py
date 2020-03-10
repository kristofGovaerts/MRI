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
import subprocess

#Globals
FUNC_DIR = "/home/brain/Desktop/TP_fMRI/EPI-rsfMRI_300kHz/rs_fMRI_MC_SC"
EPI_ana_DIR = "/home/brain/Desktop/TP_fMRI/EPI-rsfMRI_300kHz/BFC_anatomicalEPI"
HR_ANA_DIR = "/home/brain/Desktop/TP_fMRI/EPI-rsfMRI_300kHz/BFC_TR3D"
FUNC_FILELIST = "/home/brain/Desktop/TP_fMRI/EPI-rsfMRI_300kHz/rs_fMRI_MC_SC/fmri_files.txt"
EPIana_FILELIST = "/home/brain/Desktop/TP_fMRI/EPI-rsfMRI_300kHz/BFC_anatomicalEPI/EPIana_files.txt"
ANA_FILELIST = "/home/brain/Desktop/TP_fMRI/EPI-rsfMRI_300kHz/BFC_TR3D/TR3D_files.txt"
ATLAS = "/home/brain/Desktop/TP_fMRI/Atlas/BNLref_DS_ax20.hdr"

VOLUME_INDEX = 0 #which volume to register
NEUROMORPH_LOC = "C:/Program Files/NeuroMorph/Packages/MIC/Applications/Modules/Macros/KULNeuroMorph"

#default pars
THRESH=0
FTYPE='float'


#functions
def flirt(flt, atl, rdir=os.getcwd(), mode='A'):
    '''Registers [atl] into [flt] space and resamples [flt] into [atl] space.
    Inputs should be image paths.
    Can do affine (mode='A') or rigid (mode='R') registration.'''
    omat = os.path.join(rdir, os.path.basename(flt.split('.')[0] + mode + '.mat'))
    out = os.path.join(rdir, os.path.basename(flt.split('.')[0] + mode + '.nii.gz'))
    if mode == 'A': #affine    
        args = '-in %s -ref %s -out %s -omat %s' %(flt, atl, out, omat)         
    elif mode == 'R': #rigid
        args = '-in %s -ref %s -rigOnly -out %s -omat %s' %(flt, atl, out, omat)   
        
    flirt = 'flirt '+ args  
 
    process = subprocess.call(flirt, shell=True)
    return out, omat
    
def fnirt(flt, atl, rdir=os.getcwd(), aff=None, ext=''):
    '''Registers [atl] into [flt] space and resamples [flt] into [atl] space.
    Inputs should be image paths.
    Does nonrigid registration.'''
 
    out = os.path.join(rdir, os.path.basename(flt.split('.')[0] + 'NRR' + str(ext) +  '.nii')) 
    dfout = os.path.join(rdir, os.path.basename(flt.split('.')[0] + 'CPP' + str(ext) + '.nii')) 
    
    args = '--in=%s --ref=%s --cout=%s --iout=%s --subsamp=4,3,2 --miter=3,4,5 --infwhm 2,2,2 --reffwhm 2,0,0' %(flt, atl, dfout, out) 
    if aff != None:
        args += " --aff=" + aff
    fnirt = 'fnirt ' + args      
    print fnirt
    subprocess.call(fnirt, shell=True)
    return out, dfout

def combine_transforms(t1, t2, ref1, ref2, filename):
    '''Inputs:
    t1, t2: Filenames of transforms to combine (Should point to control point imgs)
    ref: Reference image used for one of the transforms
    filename: Filename for the output control point image
    Output:
    Filename of the combined transform control point image'''
    args = '-comp %s %s %s -ref %s -ref2 %s' %(t1, t2, filename, ref1, ref2)     
    tf = [reg_transform] + args.split()  
    subprocess.call(tf, shell=True)
    return filename
    
def apply_transform(t, ref, flo, filename):
    '''Inputs:
    t: Filename of a control point image or affine transformation matrix
    ref: Reference image
    flo: Floating image to which to apply the transform
    filename: Name to give the transformed floating image.
    Output:
    Filename of the transformed floating image.'''
    args = '-ref %s -flo %s -trans %s -res %s' %(ref, flo, t, filename)     
    tf = [reg_resample] + args.split()  
    print tf
    subprocess.call(tf, shell=True)
    return filename    
    
def get_deformation(t, ref, filename):
    '''Inputs:
    t: Filename of a control point image or affine transformation matrix
    ref: Reference image
    flo: Floating image to which to apply the transform
    filename: Name to give the transformed floating image.
    Output:
    Filename of the transformed floating image.'''
    args = '-def %s %s -ref %s' %(t, filename, ref)     
    tf = [reg_transform] + args.split()  
    print tf
    subprocess.call(tf, shell=True)
    return filename  
    
def nr_convert(im, out):
    args = '-in %s -out %s -float' %(im, out)     
    com = [reg_tools] + args.split()  
    subprocess.call(com, shell=True)
    return out

def batch_combine_apply_tf(ref, flo_list, tf1_list, tf2_list, ext='CPPC.nii'):
    '''Combines transforms from lists 1 and 2 and applies them to images in flo_list. 
    All input lists should be equally long. ref is the reference image used to 
    generate the transforms.'''
    ext2=ext.split('.')[0] + 'T.nii'
    tflist = []
    tfimlist = []
    assert len(flo_list) == len(tf1_list) and len(flo_list) == len(tf2_list) #inputs should be the same size!    
    
    for i, f in enumerate(flo_list):
        print "Combining transforms for image %s of %s..." %(i+1, len(flo_list))
        fn=combine_transforms(tf1_list[i], tf2_list[i], ref, f.split('.')[0] + ext)
        tflist.append(fn)
        print "Applying transform for image %s of %s..." %(i+1, len(flo_list))
        fn2=apply_transform(fn, ref, f, f.split('.')[0] + ext2)
        tfimlist.append(fn2)
    return tflist, tfimlist
    
def nrr_multiple_iterations(affine_temp, affine_list, n=5, thresh=0, ftype='int32'):
    print "ITERATION 1: Registering affine transformed images to affine template."
    if len(glob.glob('*NRR1.nii')) == len(affine_list):
        print "Registration already performed. Skipping iteration 1."
        nrr_list1 = glob.glob('*NRR1.nii')
        tf_list1 = glob.glob('*CPP1.nii')
    else:
        nrr_list1, tf_list1 = reg_batch(affine_temp, affine_list, mode='NRR', ext='1')   
        
    print "Updating template: Averaging images..."
    nrr_temp = 'nrr_template1.nii'
    save_sum_img(nrr_list1, nrr_temp, thresh=thresh, ftype=ftype)
    
    for i in range(n-1):
        print "ITERATION %s: Registering nonrigidly transformed images to template." %(i+2)
        e = '*NRR%s.nii' %(i+2)
        e2 = '*CPP%s.nii' %(i+2)
        if len(glob.glob(e)) == len(affine_list):
            print "Registration already performed. Skipping iteration %s." %(i+2)
            nrr_list2 =  glob.glob(e)
            tf_list2 = glob.glob(e2)
        else:
            nrr_list2, tf_list2 = reg_batch(nrr_temp, nrr_list1, mode='NRR', ext=(i+2))
        
        print "Updating template: Combining iterations."
        tflist, tfimlist = batch_combine_apply_tf(nrr_temp, affine_list, tf_list1, tf_list2)
        
        print "Averaging images..."
        nrr_temp='nrr_template%s.nii' %(i+2)
        save_sum_img(tfimlist, nrr_temp, thresh=thresh, ftype=ftype)
        
        tf_list1 = tflist #Keep updating your average transforms
        nrr_list1 = tfimlist
    return nrr_temp, tflist, tfimlist

    
    
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
    
def read_list(filelist, joindir = None):
    with open(filelist) as f:
        outlist = f.readlines()
    outlist = [os.path.join(joindir, l.replace('\n', '').replace('\r','')) for l in outlist]
    return outlist
    
def resample_functional(func, tf, ref, full_shape, out):
    '''Resamples functional image using a recognized transform - this is a 
    fairly slow process because it needs to resample each volume seperately
    and then compile the image again.'''
    fim = nib.load(func)
    affine = fim.get_affine()
    start_pos = -((full_shape-fim.shape[-2])/2) * affine[2,2] #Home position for expanded EPI is not the same as that for the functional img
    affine2 = np.array(affine)
    affine2[:2,3] = 0.0
    affine2[2,3]=start_pos
    rslist = []
    for t in range(fim.shape[-1]):
        fn = os.path.join(rootdir, "temp%s.nii" %t)
        fn2 = os.path.join(rootdir, "temp%sRS.nii" %t)
        d = fim.get_data()[:,:,:,0].reshape(list(fim.shape[:-1]) + [1])
        nif = nib.Nifti1Image(fim.get_data()[:,:,:,0], affine)
        nib.save(nif, fn)
        rslist.append(apply_transform(tf, ref, fn, fn2))

#Main loop
#print "Enter the scan directory. This directory should contain all scans except for the reference image, and no other scans."
#root = Tkinter.Tk()
#while True:
#    rootdir = tkFileDialog.askdirectory(initialdir="/",title='Please select a directory')
#    if os.path.isdir(rootdir) is True: #Checks if entered dir exists
#        os.chdir(rootdir)
#        root.destroy()
#        break
#    else:
#        print "Pathname invalid. Try again."
#        continue
    
rootdir = '/home/brain/Desktop/TP_fMRI/EPI-rsfMRI_300kHz/Reg_to_atlas'
os.chdir(rootdir)

funcfiles = read_list(FUNC_FILELIST, joindir=FUNC_DIR)
epifiles = read_list(EPIana_FILELIST, joindir=EPI_ana_DIR)
anafiles = read_list(ANA_FILELIST, joindir=HR_ANA_DIR)

if len(funcfiles) != len(epifiles) or len(funcfiles) != len(anafiles): #check!
    raise ValueError("Error: Input lists not of same length.")
    
epifiles = [epifiles[0]] #FOR TESTING - REMOVE

print "STEP 1: Nonrigid reg of expanded EPI image to anatomical image."
for n, f in enumerate(epifiles):
    print "File %s of %s." %(n+1, len(epifiles))
    
    epi = f
    ana = anafiles[n]
    func = funcfiles[n]
    
    print "Registering expanded EPI to anatomical image."
    aff, mat = flirt(epi, ana, rdir=rootdir, mode='A')
    nrr, cpp1 = fnirt(epi, ana, rdir=rootdir, aff=mat)
    
    print "Registering anatomical image to atlas."
    aff2, mat2 = flirt(ana, ATLAS, rdir=rootdir, mode='A')
    nrr2, cpp2 = fnirt(ana, ATLAS, rdir=rootdir, aff=mat2)