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

#Globals
N_ITERATIONS = 5 
VOLUME_INDEX = 0 #which volume to register
NEUROMORPH_LOC = "C:/Program Files/NeuroMorph/Packages/MIC/Applications/Modules/Macros/KULNeuroMorph"
#default pars
BENDING_ENERGY = 0.05
GRID_SPACING = -5 #final grid spacing - if negative, is in voxels, if positive in mm
THRESH=0
FTYPE='float'
REMOVE_BG = True #for normalized data , background value may not be 0
MODE='R'
#file extensions
tp1_ext = '*_a*FSN1.hdr'
tp2_ext = '_TR3D_2*FSN1.hdr'
out_ext = '_mid_a2.nii'
lbl_ext = '_Vent.nii'

reg_ala = os.path.join(NEUROMORPH_LOC, 'reg_aladin.exe')
reg_f3d = os.path.join(NEUROMORPH_LOC, 'reg_f3d.exe')
reg_transform = os.path.join(NEUROMORPH_LOC, 'reg_transform.exe')
reg_resample = os.path.join(NEUROMORPH_LOC, 'reg_resample.exe')
reg_jacobian = os.path.join(NEUROMORPH_LOC, 'reg_jacobian.exe')
reg_average = os.path.join(NEUROMORPH_LOC, 'reg_average.exe')

#functions
def RegAndResample(flt, atl, mode='A'):
    '''Registers [atl] into [flt] space and resamples [flt] into [atl] space.
    Inputs should be image paths.
    Can do affine (mode='A') or rigid (mode='R') registration.'''
    omat = os.path.join(os.getcwd(), flt.split('.')[0] + mode + '.txt')     
    out = os.path.join(os.getcwd(), flt.split('.')[0] + mode + '.nii') 
    if mode == 'A': #affine    
        args = '-flo %s -ref %s -noSym' %(os.path.join(os.getcwd(), flt), atl)         
    elif mode == 'R': #rigid
        args = '-flo %s -ref %s -noSym -rigOnly' %(os.path.join(os.getcwd(), flt), atl)   
        
    flirt = [reg_ala] + args.split()    
 
    process = subprocess.call(flirt, shell=True)
    os.rename('outputResult.nii', out) #NiftyReg saves as outputResult by default
    os.rename('outputAffine.txt', omat)

    return out
def otsu_filter(img, background=0):
    flat = (img.flatten())[img.flatten()>background]
    thresh = filter.threshold_otsu(flat)
    im2 = np.copy(img)
    im2[im2<=thresh]=0
    return im2
    
def reg_batch(refim, filelist, mode, zoom=1, ext='', ln=3, sx=-5, be=0.005):
    '''Performs affine or rigid registration of 'filelist' to 'refim' and 
    resamples into 'refim' space. 'mode' should be 'rigid' or 'affine'.'''
    outlist = []
    dfoutlist = []
    for i, f in enumerate(filelist):
        print "file %s of %s." %(i+1, len(filelist))

        if mode == 'NRR':
            out, dfout = RegAndResample_NRR(f, refim, ext, ln=ln, sx=sx, be=be)
            dfoutlist.append(dfout)
        else:
            out = RegAndResample(f, refim, mode=mode) #one registration step
        outlist.append(out)
    if mode == 'NRR':
        o=(outlist, dfoutlist)
    else:
        o=outlist
    return o

def combine_transforms(t1, t2, ref, filename):
    '''Inputs:
    t1, t2: Filenames of transforms to combine (Should point to control point imgs)
    ref: Reference image used for one of the transforms
    filename: Filename for the output control point image
    Output:
    Filename of the combined transform control point image'''
    args = '-comp %s %s %s -ref %s' %(t1, t2, filename, ref)     
    tf = [reg_transform] + args.split()  
    subprocess.call(tf, shell=True)
    return filename
    
def apply_transform(t, ref, flo, filename):
    '''Inputs:
    t: Filename of a control point image 
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
        fn=combine_transforms(tf2_list[i], tf1_list[i], ref, f.split('.')[0] + ext)
        tflist.append(fn)
        print "Applying transform for image %s of %s..." %(i+1, len(flo_list))
        fn2=apply_transform(fn, ref, f, f.split('.')[0] + ext2)
        tfimlist.append(fn2)
    return tflist, tfimlist
    
def save_sum_img(filelist, filename, thresh=0, ftype='int32'):
    cim = nib.load(filelist[0]).get_data()
    sum_img = np.zeros((cim.shape[0], cim.shape[1], cim.shape[2], len(filelist)))
    
    for i, f in enumerate(filelist):
        cim = np.squeeze(nib.load(f).get_data())
        sum_img[:,:,:,i] = cim  
       
    nif1 = nib.Nifti1Image(sum_img, nib.load(filelist[0]).get_affine())
    nib.save(nif1, filename[:-4] + 'stack.nii')       
     
    adata = (np.mean(sum_img, axis=-1)).astype(ftype) 
    adata[adata<thresh] = 0 #for intensity normalized images, the background is usually not 0. 
    nif = nib.Nifti1Image(adata, nib.load(filelist[0]).get_affine())
    nib.save(nif, filename)
    
def invert_transform(tf, out, template=None, fl=None):
    '''
    Inputs:
    tf: Control point image of the transform you want to invert. Can also be an affine matrix file.
    fl: Floating image used for the registration.
    template: Template image used for the registration.
    out: Output filename.
    
    Output:
    out: Output filename.
    '''
    if tf[-4:] == '.nii':
        print "Inverting nonrigid transform."
        args = '-invNrr %s %s %s -ref %s' %(tf, fl, out, template)
        reg_t = [reg_transform] + args.split()      
        subprocess.call(reg_t, shell=True)
    elif tf[-4:] == '.txt':
        print "Inverting affine transform."
        args = '-invAff %s %s ' %(tf, out)
        reg_t = [reg_transform] + args.split()      
        subprocess.call(reg_t, shell=True)
    else:
        print "Cannot identify transform type."
    return out

def batch_inverseAff(tf_list):
    in_list=[]
    for t in tf_list:
        print "File", t
        print "Computing inverse transform..."  
        out = t.replace('.txt', 'inv.txt')
        in_list.append(invert_transform(t, out))
    return in_list
    
def resave_zoomed(im, remove_bg = False, zf = 1.):
    data = nib.load(im).get_data()
    if remove_bg:
        data[data==np.min(data)]=0
    affine = nib.load(im).get_affine()
    aff2 = affine * np.diag([zf,zf,zf,1.])
    aff2[:,3] = affine[:,3]
    nif = nib.Nifti1Image(data, aff2)
    nib.save(nif, im[:-4] + 'N.nii')
    return im[:-4] + 'N.nii'
    
def aff_longitudinal(tp1_im, tp2_ext):
    tp2_list = glob.glob(tp1_im[:tp1_im.find(tp1_ext[1:3])] + tp2_ext) #find extension in filename, then use everything before the extension as root
    if len(tp2_list) > 1:
        print tp2_list
        tp2_im = tp2_list[-1]
        print "Multiple follow-up scans found. Using %s." %tp2_im
    else:
        tp2_im = tp2_list[0]
        print tp2_im
        
    tp1 = resave_zoomed(tp1_im, remove_bg = REMOVE_BG)
    tp2 = resave_zoomed(tp2_im, remove_bg = REMOVE_BG)
    tp2to1 = RegAndResample(tp2, tp1, mode=MODE)
    tp2to1_mat = tp2to1.replace('.nii', '.txt')
    
    return tp2, tp2to1, tp2to1_mat
    
def to_halfspace(tp1, tp2, tp2to1_mat, out_ext, lbl1 = None, lbl2 = None):
    forward_half = tp2to1_mat[:-4] + 'half.txt'
    backward_half = tp2to1_mat[:-4] + 'halfI.txt'
    t1m = tp1[:-4] + out_ext
    t1lbl = tp1[:-4] + '_Vent' + out_ext
    t2m = tp2[:-4] + out_ext
    t2lbl = tp2[:-4] + '_Vent' + out_ext
    print "Getting half forward and backward transforms..."
    #get half forward transform    
    args = '-half %s %s' %(tp2to1_mat, forward_half)     
    tf = [reg_transform] + args.split()  
    subprocess.call(tf, shell=True)
    
    #get half backwards transform
    args = '-invAff %s %s' %(forward_half, backward_half)     
    tf = [reg_transform] + args.split()  
    subprocess.call(tf, shell=True)
    
    print "Resampling..."
    #resample imgs
    apply_transform(forward_half, tp1, tp2, t2m)
    apply_transform(backward_half, tp2, tp1, t1m)  
    if lbl1 != None:
        apply_transform(forward_half, tp1, lbl2, t2lbl)
    if lbl2 != None:
        apply_transform(backward_half, lbl2, lbl1, t1lbl) 
    return forward_half, backward_half, t1m, t2m, t1lbl, t2lbl
    
def imgs_to_halfspace(filelist, tp2_ext):
    for i, f in enumerate(filelist):
        print "Image %s of %s." %(i+1, len(filelist))
        print "Baseline img:", f        
        try:
            tp2, tp2to1, tp2to1_mat = aff_longitudinal(f, tp2_ext)
        except IndexError:
            print "Followup image not found."
            continue
        print "Follow-up img:", tp2
        lbl1 = f[:-4] + lbl_ext
        lbl2 = tp2[:-5] + lbl_ext
        fw_half, bw_half, t1m, t2m, t1lbl, t2lbl = to_halfspace(f[:-4] + 'N.nii', tp2, tp2to1_mat, out_ext, lbl1 = lbl1, lbl2 = lbl2)
        print "Baseline output:", t1m
        print "Follow-up output:", t2m
        
        


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

filelist_tp1=glob.glob(tp1_ext)
if len(filelist_tp1) == 0:
    print "Did not find any files."
    
imgs_to_halfspace(filelist_tp1, tp2_ext)