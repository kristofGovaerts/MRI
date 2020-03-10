# -*- coding: utf-8 -*-
"""
Created on Mon May 19 09:47:12 2014

@author: brain
"""

import os
from ReadInput import *
from ImProcessing import *
from Saveoutput import *
from imageTypes import *
import nibabel as nib
import pylab as pl
import numpy as np
from dipy.segment.mask import median_otsu
from scipy.special import gamma
import dipy.reconst.dti as dti


os.chdir('/media/sf_host/data/kristof/Diffusion')
dif=DiffusionImg('FVL_25117_NC_2m_a_lf1_7_1')
dif.pdata=nib.load('FVL_25117_NC_2m_a_lf1_7_1_EC.nii.gz').get_data()

b0_mask, mask = median_otsu(dif.pdata, 3, 2)
dif.pdata=b0_mask

#for RESTORE fitting
im=nib.load('FVL_25117_NC_2m_a_lf1_7_1_EC.nii.gz').get_data()
mean_std =np.mean(np.std(im[..., :dif.nA0], -1)) #variance across b0s
n = np.sum(dif.nA0)
bias = mean_std*(1. - np.sqrt(2. / (n-1)) * (gamma(n / 2.) / gamma((n-1) / 2.)))
sigma = mean_std + bias

dti_restore = dti.TensorModel(dif.gtab,fit_method='RESTORE', sigma=sigma/2.0)
fit_restore_noisy = dti_restore.fit(dif.pdata, mask=mask)

pl.subplot(1,2,1)
pl.imshow(dif.tenfit.fa[:,:,32])
pl.subplot(1,2,2)
pl.imshow(fit_restore_noisy.fa[:,:,32])


#for CSD
from dipy.reconst.csdeconv import auto_response

response, ratio = auto_response(dif.gtab, dif.pdata, roi_center=[21,44,32], roi_radius=2, fa_thr=0.65)

from dipy.viz import fvtk

ren = fvtk.ren()

evals = response[0]

evecs = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]).T

from dipy.data import get_sphere

sphere = get_sphere('symmetric724')

from dipy.sims.voxel import single_tensor_odf

response_odf = single_tensor_odf(sphere.vertices, evals, evecs)

response_actor = fvtk.sphere_funcs(response_odf, sphere)

fvtk.add(ren, response_actor)

print('Saving illustration as csd_response.png')
fvtk.record(ren, out_path='csd_response.png', size=(200, 200))