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
from ProcessImages.ReadInput import *
from ProcessImages.ImProcessing import *
from ProcessImages.Saveoutput import *
import nibabel as nib
import pylab as pl
import numpy as np
from scipy.optimize import curve_fit
from datetime import datetime
import glob
import Tkinter, tkFileDialog
root = Tkinter.Tk()

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

l1 = "VisuAcquisitionProtocol=" #uses the Try commands underneath to load this parameter
l2 = "EffectiveEchoTime="
l3 = "VisuCoreDataSlope="

for ind, file in enumerate(filelist):
    startTime = datetime.now()
    print "\nFile %r of %r. Filename: %s" %(ind +1, len(filelist), file)
    try:
        scantype = read_line(l1, file)
        if scantype != None:
            print "Scan type: ", scantype
        else:
            print "VisuAcquisitionProtocol field empty."
            scantype = "None"
    except AttributeError:
        scantype = "Empty"
        print "VisuAcquisitionProtocol field empty."

    try:
        te = list_values(read_line(l2, file))
        print "Echo times found: ", te
    except AttributeError:
        te = []
        print "EffectiveEchoTime field empty."

    try:
        slope = list_values(read_line(l3, file))
        print "%r VisuCoreDataSlope values found." %len(slope)
    except AttributeError:
        slope = [1]
        print "VisuCoreDataSlope field empty."
    try:
        slopes=list(slope)
        for i in range(len(slope)):
            if i == 0:
                pass
            elif slope[i] == slope[i-1]:
                slopes.remove(slope[i])  #removes duplicates from list

        img = nib.load(file + '.img')
        imai = img.get_data() #initial image array which still needs the slope

        if len(slopes) == 1:
            ima = imai*slopes[0]
            print "All slope values identical. Image matrix multiplied with slope."
        elif len(slopes) > 1:
            ima = np.array(ima)  #need to know what to do if not all values are equal
            print "Slope values are not identical."

        print "Image shape: ", img.shape

        if len(te)>1 and file[-2:] == '_1':  #makes sure it is not a processed image
            print "Can do T2 relaxometry on this image. Proceeding..."
            t2img=T2Img(ima, te)
            t2map, stdevmap, si0map, soffmap, t2fig=T2Img.processimage(t2img)
            if os.path.isdir(file) is False:
                os.mkdir(file)
            np.save(os.path.join(file, file + '_t2map'), t2map)
            np.save(os.path.join(file, file + '_si0map'), si0map)
            np.save(os.path.join(file, file + '_soffmap'), soffmap)
            np.save(os.path.join(file, file + '_stdevmap'), stdevmap)
            saveslices3d(t2map, [0.0,75.0], os.path.join(file, file + "_T2map")) #saves .png of all slices in one figure
            saveslices3d(stdevmap, [0.0,10.0], os.path.join(file, file + "_Std Devs"))
            savetiffs(t2map, file) #save .tiff for each slice for further processing
            t2fig.savefig(os.path.join(file, file + "_T2fit"))
            del t2map, stdevmap, si0map, soffmap, t2img
            pl.close()
        elif 'fair' in scantype.lower() and file[-2:] == '_1': #makes sure it is not a processed image
            print "Can do FAIR-RARE ASL processing on this image. Proceeding..."
            ti=[300, 500, 700, 900, 1100, 1300, 1500, 1700, 2000, 2300, 2700, 3000, 3500, 4000]
            t1bl = 2400.0
            aslimg=ASLImg(ima, ti, t1bl)
            aslimg.processimage()
            if os.path.isdir(file) is False:
                os.mkdir(file)
            np.save(os.path.join(file, file + '_selt1'), aslimg.selt1)
            np.save(os.path.join(file, file + '_nselt1'), aslimg.nselt1)
            np.save(os.path.join(file, file + '_selstd'), aslimg.selstd)
            np.save(os.path.join(file, file + '_nselstd'), aslimg.nselstd)
            np.save(os.path.join(file, file + '_CBF'), aslimg.cbfmap)
            #saveslices3d(aslimg.cbfmap, [-100.0, 600.0], file + "_CBFmap")
            savetiffs(aslimg.cbfmap, file)
            aslimg.selt1fig.savefig(os.path.join(file, file + "_selT1fit"))
            aslimg.nselt1fig.savefig(os.path.join(file, file + "_nonselT1fit"))
            pl.close()
        else:
            print "Can't process this image, sorry."

        timecomp = datetime.now()-startTime

        print "Process took %rs to complete." %timecomp.seconds
    except NameError:
        print "NameError. No text file found or image is already processed. Skipping."
        continue
    except MemoryError:
        print "MemoryError. Scan is too large to process."
        continue
    del imai, ima, img
print "Processed all images."