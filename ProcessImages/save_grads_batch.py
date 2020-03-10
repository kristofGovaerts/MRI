####
# Copyright 2013 by Kristof Govaerts <kristof88@gmail.com>
#
# Feel free to modify any part of the code below and to contact me with
# any questions. This tool was developed at KU Leuven, Biomedical MRI Unit.
###

'''
Program will iterate over a list of files, check if they are appropriate for
T2 relaxometry analysis and process as needed.
This package is dependent on SciPy, NumPy and NiBabel. Errors can occur if
you do not have these installed.
'''

import os
from ReadInput import *
from ImProcessing import *
from Saveoutput import *
from imageTypes import *
import nibabel as nib
import pylab as pl
import numpy as np
from scipy.optimize import *
from datetime import datetime
import glob
import Tkinter, tkFileDialog
root = Tkinter.Tk()
import warnings

#Globals
REMOVEA0=2
SAVEBM=True
SAVEGRADS=True
SPINLOCKDUR = [10,30,50,70,90] #for T1rho mapping - adjust if necessary
T1VTR = [220, 4000, 350, 2000, 500, 1000] #For T1VTR mapping - adjust if necessary
MASK=True #do you want to use a brain mask? this works for many types of images but may fail - if you are missing crucial voxels then set this to False. Only change this for T2 and T1rho data.

print "Enter the scan directory."

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
print "Found %r Analyze images in directory." %len(glob.glob('*.img'))
print "Scan types found: ", list_scans()

filelist = [x.replace('.img', '') for x in glob.glob('*.img')] #cool list comprehension that gets all files

for ind, filen in enumerate(filelist):
    startTime = datetime.now()
    print "\nFile %r of %r. Filename: %s" %(ind +1, len(filelist), filen)
    try:
        im=Br2AInfo(filen)
    except AttributeError:
        print "Error. No text file found."
        continue
    except IOError:
        warnings.warn("Cannot find .hdr or .img file.")
        continue

    print "Scan type: ", im.protocol
        
    if ("dti" in im.protocol.lower() or "dki" in im.protocol.lower() or 'trace' in im.protocol.lower()) and im.name[-1]=='1':
        dtimg=DiffusionImg(filen)
        if not (len(dtimg.bvals)>=1 and len(dtimg.dwdir)>1):
            print "Diffusion image found, but insufficient diffusion dirs or b-vals."
            continue
        if os.path.isfile(dtimg.name+'_EC.nii.gz') or os.path.isfile(dtimg.name+'_EC.nii'):
            print "Found EC file. Not running new eddy current correction."
            if os.path.isfile(dtimg.name+'_EC.nii.gz'):
                dtimg.pdata=nib.load(dtimg.name+'_EC.nii.gz').get_data()
            elif os.path.isfile(dtimg.name+'_EC.nii'):
                dtimg.pdata=nib.load(dtimg.name+'_EC.nii').get_data()
            try:
                dtimg.tensorFit(bv=None, removea0=REMOVEA0, m='dti')
                dtimg.compose_image()
                dtimg.save_output(filen+'_EC_fa_md_ad_rd')
                dtimg.save_mask()
                dtimg.save_mrtrix_grads()
            except ValueError:
                print "Matrices not aligned."
        else:
            try:
                dtimg.processDiffusion(ec=True, bv=None, mode='dti', removea0=REMOVEA0)
                dtimg.compose_image(mode='dti')
                dtimg.save_output(filen+'_EC_fa_md_ad_rd')
                dtimg.save_mask()
                dtimg.save_mrtrix_grads()
            except ValueError:
                print "Matrices not aligned."
        pl.close()
    
    else:
        print "Can't process this image, sorry."

    timecomp = datetime.now()-startTime

    print "Process took %rs to complete." %timecomp.seconds
print "Processed all images."