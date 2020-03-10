# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 15:56:57 2015

@author: u0091609
"""

from dipy.align import aniso2iso
import glob
import nibabel as nib
import numpy as np




filelist=glob.glob('*_basic_Full_V1.nii')
temp_zooms=np.diag(nib.load('FVL_19657_1L_a_iE1_6_1_Brain_M4_FS.hdr').get_affine())
for f in filelist:
    im=nib.load(f)
    aff=im.get_affine()
    zooms=np.diag(aff)
    n_arr, n_aff=aniso2iso.resample(im.get_data(), aff, zooms[:3], temp_zooms[:3])
    nif=nib.Nifti1Image(n_arr, n_aff)
    nib.save(nif, f[:-4]+'HR.nii')