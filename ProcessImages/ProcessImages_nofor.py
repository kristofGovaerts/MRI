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
from ReadText import *
from ImProcessing import *
import nibabel as nib
import pylab as pl
import numpy as np
from scipy.optimize import curve_fit
from datetime import datetime
import glob

os.chdir("C:\Users\Kristof\Dropbox\NeurodegenerativeDisease\Python\Testfiles")
print "Found %r Analyze images in directory." %len(glob.glob('*.img'))
print "Scan types found: ", list_scans()

file = "FVL_21335_1R_b_hQ1_7_1"
filelist = [x.replace('.img', '') for x in glob.glob('*.img')] #cool list comprehension that gets all files

l1 = "VisuAcquisitionProtocol=" #uses the Try commands underneath to load this parameter
l2 = "EffectiveEchoTime=" 
l3 = "VisuCoreDataSlope="


startTime = datetime.now()
print file
print "File %r of %r." %(ind +1, len(filelist))
try: 
    scantype = read_line(l1, file)
    print "Scan type: ", scantype
except AttributeError:
    scantype = "none"
    print "VisuAcquisitionProtocol field empty."
    
try:
    te = list_values(read_line(l2, file))
    print "Echo times found: ", te
except AttributeError:
    te = [0]
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

    if len(te)>1:
        print "Can do T2 relaxometry on this image. Proceeding..."
        t2img=T2Img(ima, ima.shape, te)
        t2map, si0map, soffmap=T2Img.processimage(t2img)
        np.save('t2map', t2map)
        np.save('si0map', si0map)
        np.save('soffmap', soffmap)
    else: 
        print "Can't process this image, sorry."

    timecomp = datetime.now()-startTime

    print "Process took %rs to complete." %timecomp.seconds
except NameError:
    print "No text file found. Skipping."
    continue
print "Processed all images."