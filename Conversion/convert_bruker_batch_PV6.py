# -*- coding: utf-8 -*-
"""
Created on Wed Jan 03 10:26:35 2018

@author: Kristof Govaerts: kristof.govaerts@kuleuven.be / kristof88@gmail.com

Converts all scans  from all subjects to Nifti format. Outputs text file with 
parameters from method and visu_pars files.
"""

import os
import numpy as np
import re
import nibabel as nib
import Tkinter, tkFileDialog
root = Tkinter.Tk()

MODE = 'PV6'

#This is a dictionary with known data types. Important for reading the file correctly. 
#May need to add extras if you get erors. Important is that the Bruker naming convention gets assigned a Python filetype.
dtypes = { 
'_16BIT_SGN_INT': np.int16,
'_32BIT_SGN_INT': np.int32
}

line_order= ['Filename',
             'Directory', 
             'VisuCoreByteOrder',
             'VisuAcquisitionProtocol',
             'VisuCoreDim',
             'NRepetitions',
             'NAverages',
             'EchoTime',
             'RepetitionTime',
             'EffectiveEchoTime',
             'FAIR_TI',
             'NSlices',
             'DwNDiffDir',
             'DwNDiffExpEach',
             'DwAoImages',
             'DwDir',
             'DwNDiffExp',
             'DwMaxBval',
             'DwBvalEach',
             'DwEffBval',
             'DwGradVec',
             'CEST_freqs',
             'VisuCoreFrameCount',
             'VisuCoreSize',
             'VisuCoreExtent',
             'VisuCoreFrameThickness',
             'VisuCoreUnits',
             'VisuCoreResolution',
             'VisuCoreWordType',
             'VisuCoreDataOffs',
             'VisuCoreDataSlope',
             'RG',
             'VisuCoreDataMin',
             'VisuCoreDataMax',
             'RECO_map_mode',
             'GradOrientation',
             'SliceOrientation',
             'ReadOrientation',
             'VisuCoreOrientation',
             'VisuCorePosition',
             'PixelBandwidth',
             'FlipAngle',
             'Method',
             'ScanTime',
             'ReadOffset',
             'Phase1Offset',
             'Phase2Offse',
             'SliceOffset',
             'SliceGapMode',
             'SliceGap',
             'SliceDistance',
             'PVM_SpatResol',
             'RECO_transposition',
             'ACQ_dim',
             'VisuCoreResolutionMS',
             'Subject_date']

def get_folder_list(base='.'):
    '''Gets a list of all subdirectories.'''
    out = [os.path.join(base, o) for o in os.listdir(base) 
                    if os.path.isdir(os.path.join(base,o))]
    return out

def get_floatvalues(f, i):
    '''Takes one or multiple numeric values from line i in file f, and converts 
    to either a single float value or a np.1-D array of type float.'''
    try:
        floats = float(f[i][f[i].find('=') + 1 : -1])
    except ValueError: #multiple values on next line
        temp_string = str(f[i+1:])
        floats = temp_string[:temp_string.find('##') - 3]
        floats = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", floats)   #Regular expression that extracts all numbers from the string...Not sure how it works!   
        floats = np.array([float(s) for s in floats])
        if len(floats) == 1:
            floats = floats[0]
    return floats
    
def load_bruker(f, pars):
    '''Loads a Bruker 2dseq file at location f and reformats parameter file.
    Outputs an image array.'''
    
    #Step 1: Extract pars from par file
    dims = pars['VisuCoreSize']
    fc = pars['VisuCoreFrameCount']
    dtype = pars['VisuCoreWordType']
            
    if len(dims) < 3:
        dims += [1]
    size = dims + [fc] #create image shape list

    #Step 2: Load 2dseq file and reformat
    im_arr = np.fromfile(os.path.join(f, 'pdata/1/2dseq'), dtype=dtypes[dtype])
    im_arr = im_arr.reshape(size, order='F').astype('float') #Bruker files are ordered Fortran-style ('F')
    im_arr += pars['VisuCoreDataOffs']
    im_arr *= pars['VisuCoreDataSlope'] #slope correction
    
    nslices = pars['NSlices']
    nexp = fc/np.sum(nslices) #number of experiments - eg echo times, diffusion dirs * bvals etc
    try:
        if isinstance(nslices, (list, tuple, np.ndarray)):#multiple slice pckgs
            if len(set(nslices))!=1:  
                print "Amount of slices per slice package is not equal. Output image dimensions may not be correct."
            im_arr = im_arr.reshape((size[0], size[1], nexp, nslices[0] * len(nslices)), order='F') #reshape into 4Darr. order='F' important for multiple slice pckgs
        else: #1 slice pckg
            im_arr = im_arr.reshape((size[0], size[1], nexp, nslices)) #reshape into 4Darr
        im_arr = np.swapaxes(im_arr, 2, 3)#Z first, then T
    except ValueError:# for exceptions
        pass
    return im_arr
    
def get_parameters(f):
    '''Fetches parameters from visu_pars and method files. Output: parameter dict.'''
    method = os.path.join(f, 'method')
    par_dict = {"Directory" : os.path.join(os.getcwd(), f[2:])}
    with open(method) as metfile:
        met = metfile.readlines()
    vp = os.path.join(f, 'pdata/1/visu_pars') #load visu_pars
    with open(vp) as vpfile:
        vp = vpfile.readlines()
    #1. Method file
    for i, line in enumerate(met):  #not exactly the prettiest code. Writing a function for getting the values would be 50% shorter but less efficient as you'd have to go through the txt file each time. Also, the names of the parameter for the txt file doesn't always correspond with that of the method/visu_pars file      
        if '##$PVM_SpatResol' in line: #get matrix dimensions
            resol = met[i+1][:-1].split(' ')
            resol = [float(x) for x in resol]
            par_dict['VisuCoreResolution'] = resol
        if '##$PVM_SliceThick=' in line: #get matrix dimensions
            s_thick = float(met[i][met[i].find('=')+1:-1])
            par_dict['VisuCoreFrameThickness'] = s_thick
        if '##$PVM_SPackArrNSlices' in line:
            par_dict['NSlices'] = get_floatvalues(met, i).astype('int')
        if '##$CEST_FrequencyOffset' in line:
            par_dict['CEST_freqs'] = get_floatvalues(met, i)
        if '##$PVM_EchoTime' in line: #get echo time
            par_dict['EchoTime'] = get_floatvalues(met, i)
        if '##$PVM_RepetitionTime' in line: #get TR
            par_dict['RepetitionTime'] = get_floatvalues(met, i)
        if '##$MultiRepTime' in line: #get TR if there's multiple
            par_dict['RepetitionTime'] = get_floatvalues(met, i)
        if '##$PVM_NRepetitions' in line: #NReps
            par_dict['NRepetitions'] = int(get_floatvalues(met, i))      
        if '##$PVM_NAverages' in line: #NReps
            par_dict['NAverages'] = int(get_floatvalues(met, i))   
        if '##$PVM_FairTIR_Arr=' in line: #multi-TI
            par_dict['FAIR_TI'] = get_floatvalues(met, i)
        if '##$PVM_DwNDiffDir=' in line: #DTI
            par_dict['DwNDiffDir'] = get_floatvalues(met, i)
        if '##$PVM_DwNDiffExpEach=' in line:
            par_dict['DwNDiffExpEach'] = get_floatvalues(met, i)
        if '##$PVM_DwAoImages=' in line:
            par_dict['DwAoImages'] = get_floatvalues(met, i)
        if '##$PVM_DwDir=' in line:
            par_dict['DwDir'] = get_floatvalues(met, i)
        if '##$PVM_DwMaxBval=' in line:
            par_dict['DwMaxBval'] = get_floatvalues(met, i)
        if '##$PVM_DwBvalEach=' in line:
            par_dict['DwBvalEach'] = get_floatvalues(met, i)
        if '##$PVM_DwEffBval=' in line:
            par_dict['DwEffBval'] = get_floatvalues(met, i)
        if '##$PVM_DwGradVec=' in line:
            par_dict['DwGradVec'] = get_floatvalues(met, i)
        if '##$EffectiveTE=' in line:
            par_dict['EffectiveEchoTime'] = get_floatvalues(met, i)
        if '##$PVM_SPackArrSliceGap=' in line:
            par_dict['SliceGap'] = get_floatvalues(met, i)
        if '##$PVM_SPackArrSliceDistance=' in line:
            par_dict['SliceDistance'] = get_floatvalues(met, i)
        if '##$PVM_ScanTimeStr=' in line:
            par_dict['ScanTime'] = met[i+1][met[i+1].find('=') + 2 : -2]
        if '##$PVM_SPackArrReadOrient=' in line:
            par_dict['ReadOrientation'] = met[i+1][met[i+1].find('=') + 2 : -2]
            
    #2. Visu_pars file
    for i, line in enumerate(vp):
        if '##$VisuCoreSize' in line: #get matrix dimensions
            dims = vp[i+1][:-1].split(' ')
            dims = [int(x) for x in dims]
            par_dict['VisuCoreSize'] = dims
        if '##$VisuCoreFrameCount' in line: #get frame count
            fc = int(vp[i][vp[i].find('=') + 1 : -1])
            par_dict['VisuCoreFrameCount'] = fc
        if '##$VisuCoreDim=' in line: #get dimensionality
            par_dict['VisuCoreDim'] = int(vp[i][vp[i].find('=') + 1 : -1])
        if '##$VisuCoreWordType=' in line: #image bits
            dtype = vp[i][vp[i].find('=') + 1 : -1]
            par_dict['VisuCoreWordType'] = dtype
        if '##$VisuCoreDataOffs' in line: #offsets
            par_dict['VisuCoreDataOffs'] = get_floatvalues(vp, i)
        if '##$VisuCoreDataSlope' in line: #slopes
            par_dict['VisuCoreDataSlope'] = get_floatvalues(vp, i)
        if '##$VisuCoreByteOrder=' in line:
            par_dict['VisuCoreByteOrder']=vp[i][vp[i].find('=') + 1 : -1]
        if '##$VisuAcquisitionProtocol' in line:
            par_dict['VisuAcquisitionProtocol']=vp[i+1][vp[i+1].find('=') + 2 : -2]
        if '##$VisuAcqFlipAngle=' in line:
            par_dict['FlipAngle'] = get_floatvalues(vp, i)
        if '##$VisuCorePosition=' in line:
            par_dict['VisuCorePosition'] = get_floatvalues(vp, i)
        if '##$VisuCoreOrientation=' in line:
            par_dict['VisuCoreOrientation'] = get_floatvalues(vp, i)
        if '##$VisuCoreExtent=' in line:
            par_dict['VisuCoreExtent'] = get_floatvalues(vp, i)
        if '##$VisuCoreUnits=' in line:
            par_dict['VisuCoreUnits'] = vp[i+1][vp[i+1].find('=') + 1 : -1]
        if '##$VisuSubjectBirthDate=' in line:
            par_dict['Subject_date'] = vp[i+1][vp[i+1].find('=') + 2 : -2]
    return par_dict

def save_parfile(pardict, fn):
    '''Saves parameter dict into text file at filename fn.'''
    txt=''
    for l in line_order:
        try:
            item = pardict[l]
            if isinstance(item, (list, tuple, np.ndarray)):#formatting
                item = ' '.join(str(e) for e in item)
            else:
                item=str(item)
            txt += l + '=[' + item + ']\n'
        except KeyError:
            pass
    with open(fn, 'w') as f:
        f.write(txt)
            
def make_nifti(f, pars):
    '''Makes a Nifti file out of a Bruker folder and associated parameter dict (see also get_parameters()).'''
    arr = load_bruker(f, pars)
    if len(pars['VisuCoreResolution']) == 2:
        affine = np.diag(pars['VisuCoreResolution'] + [pars['VisuCoreFrameThickness'], 1.0])
    elif len(pars['VisuCoreResolution']) == 3:
        affine = np.diag(pars['VisuCoreResolution'] + [1.0])
    else: 
        print "Dimensionality not understood."
    nif = nib.Nifti1Image(arr, affine)
    return nif
    
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

folders = get_folder_list()
print 'Found %s folders.' %len(folders)

for i, f in enumerate(folders):
    scans = get_folder_list(f)
    scans.sort()
    scans = scans[:-1]#remove last folder: AdjResult 
    print 'folder %s of %s.' %(i+1, len(folders))
    for j, s in enumerate(scans):
        if MODE == 'PV6':
            fn = f[18:] + '_' + os.path.split(s)[-1] 
        elif MODE == 'PV5':
            fn = f[2:].replace('.','_') + '_' + os.path.split(s)[-1] 
        print 'scan %s of %s. Filename: %s' %(j+1, len(scans), fn)
        try:
            pars = get_parameters(s)
            pars['Filename'] = fn
            nif = make_nifti(s, pars)
        except IOError as e:
            print "Missing essential file."
            print e
            continue
        nib.save(nif, fn + '.nii') #index 18 removes numbers and stuff from dir. second param is scan number
        save_parfile(pars, fn + '.txt')