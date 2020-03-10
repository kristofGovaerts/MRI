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

#Globals
N_ITERATIONS = 5 
VOLUME_INDEX = 0 #which volume to register
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
    
def RegAndResample_NRR(flt, atl, ext='', ln=3, sx=-5, be=0.005):
    '''Registers [atl] into [flt] space and resamples [flt] into [atl] space.
    Inputs should be image paths.
    Does nonrigid registration.'''
 
    out = os.path.join(os.getcwd(), flt.split('.')[0] + 'NRR' + str(ext) +  '.nii') 
    dfout = os.path.join(os.getcwd(), flt.split('.')[0] + 'CPP' + str(ext) + '.nii') 
    
    args = '-flo %s -ref %s -ln %s -sx %s -be %s' %(os.path.join(os.getcwd(), flt), atl, str(ln), str(sx), str(be))     
    fnirt = [reg_f3d] + args.split()  
    
    subprocess.call(fnirt, shell=True)
    os.rename('outputResult.nii', out) #NiftyReg saves as outputResult by default
    os.rename('outputCPP.nii', dfout)
    return out, dfout
    
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
    
def nrr_multiple_iterations(affine_temp, affine_list, n=5, thresh=7000, ftype='int32'):
    print "ITERATION 1: Registering affine transformed images to affine template."
    if len(glob.glob('*NRR1.nii')) == len(affine_list):
        print "Registration already performed. Skipping iteration 1."
        nrr_list1 = glob.glob('*NRR1.nii')
        tf_list1 = glob.glob('*CPP1.nii')
    else:
        nrr_list1, tf_list1 = reg_batch(affine_temp, affine_list, mode='NRR', ext='1', sx=-3, be=0.005)   #5-voxel grid spacing
        
    print "Updating template: Averaging images..."
    nrr_temp = 'nrr_template1.nii'
    save_sum_img(nrr_list1, nrr_temp, thresh=thresh, ftype=ftype)
    
    for i in range(n-1):
        e = '*NRR%s.nii' %(i+2)
        e2 = '*CPP%s.nii' %(i+2)
        if i == 0 and n>2:
            sx = -1 #3-voxel grid spacing on first iteration
            be = 0.001
        else:
            sx = -1 #1-voxel grid spacing on all other iterations
            be = 0.001
        print "ITERATION %s: Registering nonrigidly transformed images to template. Using grid spacing of %s, 3 levels." %(i+2, sx)
        if len(glob.glob(e)) == len(affine_list):
            print "Registration already performed. Skipping iteration %s." %(i+2)
            nrr_list2 =  glob.glob(e)
            tf_list2 = glob.glob(e2)
        else:
            nrr_list2, tf_list2 = reg_batch(nrr_temp, nrr_list1, mode='NRR', ext=(i+2), ln=3, sx=sx, be=be)
        
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
    nib.save(nif1, filename[:-4] + 'stack.nii')       
     
    adata = (np.mean(sum_img, axis=-1)).astype(ftype) 
    adata[adata<thresh] = 0 #for intensity normalized images, the background is usually not 0. 
    nif = nib.Nifti1Image(adata, nib.load(filelist[0]).get_affine())
    nib.save(nif, filename)
    
def invert_transform(tf, fl, template, out):
    '''
    Inputs:
    tf: Control point image of the transform you want to invert. Can also be an affine matrix file.
    fl: Floating image used for the registration.
    template: Template image used for the registration.
    out: Output filename.
    
    Output:
    out: Output filename.
    '''
    args = '-invNrr %s %s %s -ref %s' %(tf, fl, out, template)
    reg_t = [reg_transform] + args.split()      
    subprocess.call(reg_t, shell=True)
    return out
    
def get_jacobian(df, template, out):
    '''
    Inputs:
    df: Deformation field or control point image.
    template: Template image.
    out: Output filename.
    
    Output:
    out: Output filename.
    '''
    args = '-cpp %s -target %s -jac %s' %(df, template, out)
    reg_j = [reg_jacobian] + args.split()  
    subprocess.call(reg_j, shell=True)
    return out
    
    
def batch_inverseJac(tf_list):
    temp = 'nrr_template4.nii'
    in_list=[]
    jac_list=[]
    for t in tf_list:
        f = t.replace('CPPC.nii', '.nii') #floating image is necessary
        out_cpp = t[:-4] + "Inv.nii" 
        out_jac = t[:-4] + "InvJ.nii" #output deformation field
        print "File", t
        print "Computing inverse transform..."        
        in_list.append(invert_transform(t, f, temp, out_cpp))
        
        print "Computing jacobian..."
        jac_list.append(get_jacobian(out_cpp, f, out_jac))
        
    return in_list, jac_list
        


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

#rootdir = '/media/sf_host/data/NRR/Raw_data/test_FSL'
#os.chdir(rootdir)
#refim = '/media/sf_host/data/NRR/Raw_data/Bl6_Ref_ER.hdr'
#maskim = '/media/sf_host/data/NRR/Raw_data/Bl6_ref_ER_BNL_lbl_BM.hdr'

filelist=glob.glob('*.hdr')
if len(filelist) == 0:
    filelist=glob.glob('*.nii')
#print "STEP 1: Performing rigid registration of all images to reference image."
#rigid_list = reg_batch(refim, filelist, mode='R')
#print "STEP 2: Averaging registered images to create rigid template."
#rigid_temp = 'rigid_template.nii'
#save_sum_img(rigid_list, rigid_temp, thresh=THRESH, ftype=FTYPE)
print "STEP 3: Performing affine registration of all images to reference image."
affine_list = reg_batch(refim, filelist, mode='A')
print "STEP 4: Averaging registered images to create affine template."
affine_temp = 'affine_template.nii'
save_sum_img(affine_list, affine_temp, thresh=THRESH, ftype=FTYPE)
print "STEP 5: Iteratively creating template using nonrigid registration."
nrr_multiple_iterations(affine_temp, affine_list, n=N_ITERATIONS, thresh=THRESH, ftype=FTYPE)
print "STEP 6: Computing inverse transform."
transforms_list=glob.glob('*CPPC.nii')
batch_inverseJac(transforms_list)