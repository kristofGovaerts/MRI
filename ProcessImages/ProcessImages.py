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
#T1VTR = [220, 4000, 350, 2000, 500, 1000] #For T1VTR mapping - adjust if necessary
T1VTR = [270.5945, 6000, 429.323, 2522.966, 1724.751, 1249.591, 910.113, 645.798]
TI = [300, 500, 700, 900, 1100, 1300, 1500, 1700, 2000, 2300, 2700, 3000, 3500, 4000] #14TI
#TI = [50, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 2000, 2300, 2700, 3000, 3500, 4000, 8000] #16TI
#TI = [300, 500, 700, 900, 1100, 1300, 1500, 1700, 2000, 3000] #10TI_TS
MASK=True #do you want to use a brain mask? this works for many types of images but may fail - if you are missing crucial voxels then set this to False. Only change this for T2 and T1rho data.
EC=True
NPACKS=1

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

#    if len(im.te)>2 and filen[-2:] == '_1' and os.path.isdir(filen) is False:  #makes sure it is not a processed image, either by Paravision or by this script
    if len(im.te)>2 and filen[-2:] == '_1': 
        if os.path.exists(filen):
            print "Processed file directory already exists. Skipping..."
            continue
        print "Can do T2 relaxometry on this image. Proceeding..."
        t2img=T2Img(filen, nlm_denoise=False)
        t2img.npacks = NPACKS
        if t2img.npacks != 1:
            newshape = [t2img.shape[0], t2img.shape[1], t2img.shape[2] * t2img.npacks, t2img.shape[3] / t2img.npacks]
            t2img.pdata = np.reshape(t2img.pdata, newshape)
            t2img.shape = newshape
        if os.path.isdir(filen) is False:
            os.mkdir(filen)
#        pnif = nib.Nifti1Image(t2img.pdata, t2img.affine)
#        nib.save(pnif, os.path.join(filen, filen + '_denoise.nii'))
        print "Echo times found: ", t2img.te
        t2img.processImage(mask=MASK)
        t2img.processImage_basic(mask=MASK)
        np.save(os.path.join(filen, filen + '_t2map'), t2img.t2map)
        np.save(os.path.join(filen, filen + '_si0map'), t2img.si0map)
        np.save(os.path.join(filen, filen + '_soffmap'), t2img.soffmap)
        np.save(os.path.join(filen, filen + '_stdevmap'), t2img.stdevmap)
        
        np.save(os.path.join(filen, filen + '_t2map'), t2img.t2map_basic)
        np.save(os.path.join(filen, filen + '_si0map'), t2img.si0map_basic)
        np.save(os.path.join(filen, filen + '_stdevmap'), t2img.stdevmap_basic)
        saveslices3d(t2img.t2map, [0.0,75.0], os.path.join(filen, filen + "_T2map")) #saves .png of all slices in one figure
        saveslices3d(t2img.stdevmap, [0.0,10.0], os.path.join(filen, filen + "_Std Devs"))
        t2img.t2fig.savefig(os.path.join(filen, filen + "_T2fit"))
        t2img.t2fig_basic.savefig(os.path.join(filen, filen + "_T2fit_basic"))
        t2img.save_output(os.path.join(filen, filen + '_Full'))
        t2img.save_output_basic(os.path.join(filen, filen + '_basic_Full'))
        t2nif=nib.Nifti1Image(t2img.get_te1(mask=MASK), t2img.affine)
        nib.save(t2nif, os.path.join(filen,filen+'_te1.nii'))
        pl.close()
        
    elif 'vtr' in im.protocol.lower() and filen[-2:] == '_1': 
        print "Can do T1_VTR relaxometry on this image. Proceeding..."
        t2img=T1Img(filen)
        o_te = list(t2img.te)
        t2img.te = T1VTR #Conversion tool does not store variable TRs - Input manually. Easier to just swap the TE with TR
        print "Repetition times found: ", t2img.te
        print "Correcting for slope..."
        t2img.correctSlope()
        t2img.processImage(mask=MASK)
        t2img.processImage_basic(mask=MASK)
        if os.path.isdir(filen) is False:
            os.mkdir(filen)
        np.save(os.path.join(filen, filen + '_t2map'), t2img.t2map)
        np.save(os.path.join(filen, filen + '_si0map'), t2img.si0map)
        np.save(os.path.join(filen, filen + '_soffmap'), t2img.soffmap)
        np.save(os.path.join(filen, filen + '_stdevmap'), t2img.stdevmap)
        
        np.save(os.path.join(filen, filen + '_t2map'), t2img.t2map_basic)
        np.save(os.path.join(filen, filen + '_si0map'), t2img.si0map_basic)
        np.save(os.path.join(filen, filen + '_stdevmap'), t2img.stdevmap_basic)
        #saveslices3d(t2img.t2map, [0.0,75.0], os.path.join(filen, filen + "_T2map")) #saves .png of all slices in one figure
        #saveslices3d(t2img.stdevmap, [0.0,10.0], os.path.join(filen, filen + "_Std Devs"))
        t2nif=nib.Nifti1Image(t2img.t2map,t2img.affine)
        nib.save(t2nif, os.path.join(filen,filen+'_T2Map.nii'))
        t2img.t2fig.savefig(os.path.join(filen, filen + "_T2fit"))
        t2img.t2fig_basic.savefig(os.path.join(filen, filen + "_T2fit_basic"))
        t2img.save_output(os.path.join(filen, filen + '_Full'))
        t2img.save_output_basic(os.path.join(filen, filen + '_basic_Full'))
        t2img.te = o_te #put original TE back just in case - will probably try to do T2 fitting if it thinks there are more than 1 TE
        pl.close()
        
    elif 'fair' in im.protocol.lower() and filen[-2:] == '_1' and os.path.isdir(filen) is False: #makes sure it is not a processed image
        print "Can do FAIR-RARE ASL processing on this image. Proceeding..."
        aslimg=ASLImg(filen)
        print "Correcting for slope..."
        aslimg.correctSlope()
        aslimg.ti=TI
#        aslimg.ti=[100, 500, 900, 1300, 1800, 2500, 4000]
        aslimg.t1bl = 2400.0
        try:
            aslimg.processimage()            
            if os.path.isdir(filen) is False:
                os.mkdir(filen)
            np.save(os.path.join(filen, filen + '_selt1'), aslimg.selt1)
            np.save(os.path.join(filen, filen + '_nselt1'), aslimg.nselt1)
            np.save(os.path.join(filen, filen + '_selstd'), aslimg.selstd)
            np.save(os.path.join(filen, filen + '_nselstd'), aslimg.nselstd)
            np.save(os.path.join(filen, filen + '_CBF'), aslimg.cbfmap)
            #saveslices3d(aslimg.cbfmap, [-100.0, 600.0], file + "_CBFmap")
            aslnif=nib.Nifti1Image(aslimg.cbfmap, aslimg.affine, header=aslimg.img.get_header())
            nib.save(aslnif,os.path.join(filen, filen+'_CBFmap.nii'))
            aslimg.selt1fig.savefig(os.path.join(filen, filen + "_selT1fit"))
            aslimg.nselt1fig.savefig(os.path.join(filen, filen + "_nonselT1fit"))
            aslimg.save_output(os.path.join(filen, filen + '_Full'))
            pl.close()
        except ValueError:
            print "Image shape incorrect: ", aslimg.shape, "Not compatible with amount of inversion times provided."
        
    elif ("dti" in im.protocol.lower() or "dwi" in im.protocol.lower() or "dki" in im.protocol.lower() or 'trace' in im.protocol.lower()) and im.name[-1]=='1':
        dtimg=DiffusionImg(filen)
        if not ((len(dtimg.bvals)-dtimg.nA0)/len(dtimg.dwdir)>=1 and len(dtimg.dwdir)>5):
            if EC:    
                from PythonDiffusion import eddyCorrection
                eddyCorrection(dtimg, dtimg.name+'_EC', protocol='eddy_correct')
            if (len(dtimg.bvals)-dtimg.nA0)/len(dtimg.dwdir) == 1:
                print "Single b-value found, but not enough directions to fit a tensor. Calculating trace image."
                dtimg.trace_1b()
                nif = nib.Nifti1Image(dtimg.trace, dtimg.affine)
                nib.save(nif, filen + '_trace.nii')
            elif (len(dtimg.bvals)-dtimg.nA0)/len(dtimg.dwdir) > 1:
                print "Multiple b-values found, but not enough directions to fit a tensor. Fitting trace image."
                dtimg.trace_multib()
                nif = nib.Nifti1Image(dtimg.trace, dtimg.affine)
                nib.save(nif, filen + '_trace.nii')

        elif os.path.isfile(dtimg.name+'_EC.nii.gz') or os.path.isfile(dtimg.name+'_EC.nii'):
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
                dtimg.processDiffusion(ec=EC, bv=None, mode='dti', removea0=REMOVEA0, mask=MASK)
                dtimg.compose_image(mode='dti')
                dtimg.save_output(filen+'_EC_fa_md_ad_rd')
                dtimg.save_mask()
                dtimg.save_mrtrix_grads()
            except ValueError:
                print "Matrices not aligned."
        pl.close()
        
    elif (('t1rho' in im.protocol.lower() or '(modified)' in im.protocol.lower() or ''==im.protocol.lower())) and 'vtr' not in im.protocol.lower(): #t1rho scans - names are kind of weird, add 'or' clauses for your protocol if necessary 
        print "Can do T1rho relaxometry on this image. Proceeding..."
        t2img=T2Img(filen, repetition_avg = False)
        if t2img.shape[-1]==1:
            print 'Not enough frames. Continuing to next img. Need a 4D image.'
            continue
        t2img.te=SPINLOCKDUR
        print "Echo times found: ", t2img.te
        print "Correcting for slope..."
        if len(list(set(t2img.slopes))) == 1:
            t2img.pdata *= list(set(t2img.slopes))[0]
        else:
            print "Can't correct slope. Moving to next img."
            continue
        t2img.processImage(mask=MASK) #turn masking off - better if you're not working with brains
        t2img.processImage_basic(mask=MASK)
        if os.path.isdir(filen) is False:
            os.mkdir(filen)
        np.save(os.path.join(filen, filen + '_t2map'), t2img.t2map)
        np.save(os.path.join(filen, filen + '_si0map'), t2img.si0map)
        np.save(os.path.join(filen, filen + '_soffmap'), t2img.soffmap)
        np.save(os.path.join(filen, filen + '_stdevmap'), t2img.stdevmap)
        
        np.save(os.path.join(filen, filen + '_t2map'), t2img.t2map_basic)
        np.save(os.path.join(filen, filen + '_si0map'), t2img.si0map_basic)
        np.save(os.path.join(filen, filen + '_stdevmap'), t2img.stdevmap_basic)
        #saveslices3d(t2img.t2map, [0.0,75.0], os.path.join(filen, filen + "_T2map")) #saves .png of all slices in one figure
        #saveslices3d(t2img.stdevmap, [0.0,10.0], os.path.join(filen, filen + "_Std Devs"))
        t2nif=nib.Nifti1Image(t2img.t2map,t2img.affine)
        nib.save(t2nif, os.path.join(filen,filen+'_T2Map.nii'))
        t2img.t2fig.savefig(os.path.join(filen, filen + "_T2fit"))
        t2img.t2fig_basic.savefig(os.path.join(filen, filen + "_T2fit_basic"))
        t2img.save_output(os.path.join(filen, filen + '_Full'))
        t2img.save_output_basic(os.path.join(filen, filen + '_basic_Full'))
        pl.close()    
    
    
    else:
        print "Can't process this image, sorry."

    timecomp = datetime.now()-startTime

    print "Process took %rs to complete." %timecomp.seconds
print "Processed all images."