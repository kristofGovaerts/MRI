# -*- coding: utf-8 -*-
"""
Created on Thu May 21 11:45:03 2015

@author: u0091609
"""

import os
import nibabel as nib
import numpy as np
import glob
import Tkinter, tkFileDialog
import subprocess

#Globals
JAC_EXT='_FSN1MACPPCInvJ.nii' #file extension for the jacobian maps for your image
AFF_EXT = '_FSN1MA.txt'  #file extension for the affine matrices for your image
LABEL_IM = 'BNL_avgA_inv_DSc.hdr' #Label image with labels in refspace
UPPER_BOUND = 100

LAB_NAMES = ["HIP", "EC", "CP", "ACO", "Glob.Pall", "IC", "Thal", "CER", 
"Sup.Coll", "VENT", "HYPO", "INF.COL", "Cent.GRAY", "CORT", "AMYG", "OB", 
"Brainstem", "MIDBRAIN", "BF_sept", "FIM"]

reg_transform = os.path.join(NEUROMORPH_LOC, 'reg_transform.exe')

def read_affine(aff_filename):
    out_mat = np.zeros((4,4))
    with open(aff_filename) as f:
        for i, line in enumerate(f.readlines()):
            out_mat[:,i] = np.array(line.split())
    return out_mat

def quantify_labels(pim, lim, vvol):
    '''Quantifies the sum of the values in each label and multiplies by voxel volume - Returns label volumes in float space.'''
    hist = np.histogram(lim, bins=np.max(lim)) #Check how many labels there are
    values = hist[1][1:]
    print "Found %s labels." %len(values)
    lab_sums=[]
    lab_stdevs=[]
    lab_datapoints=[]
    cmask = np.zeros(lim.shape)
    cmask[pim[:,:,:]>0] = 1
    cmask[pim[:,:,:]>UPPER_BOUND] = 0    
    pim[np.isnan(pim)] = 0    
    
    pim[:,:,:] *= cmask 
    for v in values:
        print "label", v
        par_sums = np.sum((vvol * pim[:,:,:])[lim==v])#multiply by voxel volume to get real volumes
        par_stdevs = np.std((pim[:,:,:])[lim==v])     
        par_datapoints = (pim[:,:,:])[lim==v]               
        lab_sums.append(par_sums)
        lab_stdevs.append(par_stdevs)
        lab_datapoints.append(par_datapoints)
    return lab_sums, lab_stdevs, lab_datapoints
    
def give_histograms(datapoints, legend, pnames):
    fig=pl.figure(figsize=(40, 30), dpi=80)
    l=len(datapoints)
    npars=len(datapoints[0])
    x=np.sqrt(npars)
    for p in range(npars):
        pl.subplot(np.floor(x)+1, np.ceil(x), p+1)
        pl.title(pnames[p])
        for i in range(l):
            h=np.histogram(datapoints[i][p], 50)
            pl.plot(h[1][1:], h[0])
        pl.legend(legend)
    return fig        

def to_csv_line(string):
    substr=['[', ']', '(' , ')', 'memmap']
    o=str(string)
    for s in substr:
        o=o.replace(s, '')
    o=o.replace(',', ';')
    o=o.replace(".",",")
    return o

#Main loop
print "Enter the scan directory. This directory should contain all jacobian deformation maps and affine matrices."
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

filelist = glob.glob('*' + JAC_EXT)
lab_im = nib.load(LABEL_IM).get_data()

hdr = ';'.join([l + ' abs_vol(mm3);' + l + ' norm_vol(mm3);' + l + ' stdev' for l in LAB_NAMES])
txt = 'filename;' + hdr + '\n'

for f in filelist:
    jac = nib.load(f)
    aff_tf = f.replace(JAC_EXT, AFF_EXT)
    aff_mat = read_affine(aff_tf)
    aff_scale = np.abs(np.prod(np.diag(aff_mat)[:3]))
    with open(aff_tf) as x:
        aff_mat = x.readlines()
    jac_aff = jac.get_affine()
    zooms = np.diag(jac_aff)[:3]
    vox_vol = np.abs(np.prod(zooms))
    
    nvols, lstds, ldps = quantify_labels(jac.get_data(), lab_im, vox_vol)
    lvols = [l*aff_scale for l in nvols]
    
    zipped = zip(lvols, nvols, lstds)
    line = f + ';'
    for t in zipped:
        for i in t:
            line += str(i) + ';'
    txt += line + '\n'
    
with open('labels_quant.csv', 'w') as f:
    f.write(txt) 
    