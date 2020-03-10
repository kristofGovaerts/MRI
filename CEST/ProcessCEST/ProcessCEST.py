# -*- coding: utf-8 -*-
"""
Created on Sun Oct 06 12:37:36 2013

@author: Gebruiker
"""

from ProcessCEST.CESTimports import *
from ProcessImages.ReadInput import *
import math

import scipy


fieldmap=Bruker2AnalyzeImg('TD_GlutCEST_t02_lZ1_9_1')
fieldmap.data=scipy.swapaxes(fieldmap.data,1,2) #fieldmap array is ordered differently. This puts it in the same order as e.g. the voxel size

filelist = [x.replace('.img', '') for x in glob.glob('*.img')] #cool list comprehension that gets all files
cestimg=Bruker2AnalyzeImg('TD_GlutCEST_t02_lZ1_15_1')
cestarray=scansToArray(filelist,start=15)
cestimg.data=cestarray

fieldmap.resampleImage(cestimg) #resamples fieldmap to match CEST image resolution
fieldmap.data=fillEdges(fieldmap,cestimg)
fieldmapslice=getCorrespondingSlice(fieldmap,cestimg)

offsimg1, pseudofieldmap1, indexerrors1=correctB0(cestimg.data,ranges,1200,25,B0=B0)
offsimg2, pseudofieldmap2, indexerrors2=correctB0(cestimg.data,ranges,250,25,B0=B0)
offsimg3, pseudofieldmap3, indexerrors3=correctB0(cestimg.data,ranges,750,25,B0=B0)
offsimg4, pseudofieldmap4, indexerrors4=correctB0(cestimg.data,ranges,1500,25,B0=B0)

percentDiff1200=computePercentDiff(offsimg1[:,:,0,0],offsimg[:,:,0,1])
percentDiff250=computePercentDiff(offsimg2[:,:,0,0],offsimg[:,:,0,1])
percentDiff750=computePercentDiff(offsimg3[:,:,0,0],offsimg[:,:,0,1])
percentDiff1500=computePercentDiff(offsimg4[:,:,0,0],offsimg[:,:,0,1])

absDiff1200=computeDifference(offsimg1[:,:,0,0],offsimg[:,:,0,1])
absDiff250=computeDifference(offsimg2[:,:,0,0],offsimg[:,:,0,1])
absDiff750=computeDifference(offsimg3[:,:,0,0],offsimg[:,:,0,1])
absDiff1500=computeDifference(offsimg4[:,:,0,0],offsimg[:,:,0,1])

pl.subplot(2,2,1)
pl.imshow(percentDiff250)
pl.title('Percent difference\n250Hz')
pl.clim(0,100)
pl.colorbar()
pl.axis('off')
pl.subplot(2,2,2)
pl.imshow(percentDiff750)
pl.title('Percent difference\n750Hz')
pl.clim(0,100)
pl.colorbar()
pl.axis('off')
pl.subplot(2,2,3)
pl.imshow(percentDiff1200)
pl.title('Percent difference\n1200Hz')
pl.clim(0,100)
pl.colorbar()
pl.axis('off')
pl.subplot(2,2,4)
pl.imshow(percentDiff1500)
pl.title('Percent difference\n1500Hz')
pl.clim(0,100)
pl.colorbar()
pl.axis('off')

pl.subplot(2,2,1)
pl.imshow(absDiff250)
pl.title('Absolute difference\n250Hz')
pl.clim(0,5e7)
pl.colorbar()
pl.axis('off')
pl.subplot(2,2,2)
pl.imshow(absDiff750)
pl.title('Absolute difference\n750Hz')
pl.clim(0,5e7)
pl.colorbar()
pl.axis('off')
pl.subplot(2,2,3)
pl.imshow(absDiff1200)
pl.title('Absolute difference\n1200Hz')
pl.clim(0,5e7)
pl.colorbar()
pl.axis('off')
pl.subplot(2,2,4)
pl.imshow(absDiff1500)
pl.title('Absolute difference\n1500Hz')
pl.clim(0,5e7)
pl.colorbar()
pl.axis('off')

pl.subplot(2,2,1)
pl.imshow(cestimg.data[:,:,0,0])
pl.title('CEST image: pre-correction.\n2000Hz')
pl.colorbar()
pl.axis('off')
pl.subplot(2,2,2)
pl.imshow(B0[:,:,0,0])
pl.title('Fieldmap')
pl.colorbar()
pl.axis('off')
pl.subplot(2,2,3)
pl.imshow(offsimg1[:,:,0,0])
pl.title('B0-corrected image\n1200Hz')
pl.clim(0,1.6e8)
pl.colorbar()
pl.axis('off')
pl.subplot(2,2,4)
pl.imshow(offsimg1[:,:,0,1])
pl.title('B0-corrected image\n-1200Hz')
pl.clim(0,1.6e8)
pl.colorbar()
pl.axis('off')