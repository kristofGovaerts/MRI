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
print "Found %r Nifti images in directory." %len(glob.glob('*.nii'))
print "Scan types found: ", list_scans()

filelist = [x.replace('.nii', '') for x in glob.glob('*.nii')] #cool list comprehension that gets all files

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

    if isinstance(im.te, (list, tuple, np.ndarray)) and len(im.te)>2: 
        if os.path.exists(filen):
            print "Processed file directory already exists. Skipping..."
            continue
        print "Can do T2 relaxometry on this image. Proceeding..."
        t2img=T2Img(filen, nlm_denoise=False, slopecorr=False, reshape=False)
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
        
    elif isinstance(im.tr, (list, tuple, np.ndarray)) and len(im.tr)>2: 
        print "Can do T1_VTR relaxometry on this image. Proceeding..."
        t2img=T1Img(filen)
        o_te = list(t2img.te)
        t2img.te = t2img.tr
        print "Repetition times found: ", t2img.te
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
        
    elif 'fair' in im.protocol.lower():
        print "Can do FAIR-RARE ASL processing on this image. Proceeding..."
        aslimg=ASLImg(filen)
        aslimg.ti=np.squeeze(np.array([list_values(read_line("FAIR_TI=",filen))]))
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