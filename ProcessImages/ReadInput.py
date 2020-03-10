#Copyright  Copyright 2013 by Kristof Govaerts <kristof88@gmail.com>

"""Modules that read text files produced
after batch conversion to Analyze.
"""
import glob
import numpy as np
import nibabel as nib


def read_line(str, file):
    '''Reads file.txt and looks for str.
    If str is in a line, returns rest of line.'''
    try:
        f = open(file + '.txt')
        for l in f:
            if  str in l:
                value = l.replace(str, '')
                return value[:-1] #to remove /n char
        f.close()
    except IOError:
        pass

def list_values(str):
    '''Puts values between [] in list< '''
    li=str.replace('[','').replace(']','').split()
    return [float(x) for x in li]

def list_scans():
    txts = glob.glob('*.txt') #list of all txt files
    scanlist = [] #list of all scans, to be appended
    for txt in txts: #takes everything past = on line 6 and appends scanlist if not yet present
        line = read_line("VisuAcquisitionProtocol=", txt[:-4])
        if line in scanlist:
            pass
        else:
            scanlist.append(line)
    return scanlist
    
def sort_scans(l):
    list2=l[:]
    for i in range(100):
        s='_%r_' %i
        if i<10:
            s2='_00%r_' %i
            list2=[f.replace(s,s2) for f in list2]
        else:
            s2='_0%r_' %i
            list2=[f.replace(s,s2) for f in list2]
    list2.sort()
    for i in range(100):
        if i<10:
            s='_00%r_' %i
            s2='_%r_' %i
            list2=[f.replace(s,s2) for f in list2]
        else:
            s='_0%r_' %i
            s2='_%r_' %i
            list2=[f.replace(s,s2) for f in list2]
    return list2


def vector_array(list, m):
    '''Reshapes list of n*m values into an n x m array of vectors.'''
    va=np.array(list).reshape(len(list)/m, m)
    return va
    
def correctSlope(img, filename, pr=True):
    st="VisuCoreDataSlope="
    cimg=np.array(img)
    slopes=list_values(read_line(st, filename))
    if len(slopes)==img.shape[-1]:
        if pr:
            print "Correcting slopes..."
        for i, slope in enumerate(slopes):
            cimg[:,:,:,i]*=slope
    return cimg

def scansToArray(filelist, start=0, end=None, correct=True):
    '''Returns an array concatenating a series of 3D arrays across the
    t-dimension.

    Inputs:
        filelist = list of files, no extension
        start = where to start concatenating, optional (necessary if matrix size not equal for all images in your filelist)
        end=where to stop'''
    filelist = [x.replace('.img', '') for x in glob.glob('*.img')]
    sfl=sort_scans(filelist)
    if end is None:
        scans=sfl[start:]
    else:
        scans=sfl[start:end]
    data1=nib.load(scans[0]+'.img').get_data()
    if correct==True:
        data1=correctSlope(data1, scans[0])
    array=np.zeros([data1.shape[0], data1.shape[1],1, len(scans)])
    array[:,:,:,0]=data1[:,:,:,0]
    for i, scan in enumerate(scans[1:]):
        data=nib.load(scan+'.img').get_data()
        if correct==True:
            data=correctSlope(data, scan, pr=False)
        array[:,:,:,i+1]=data[:,:,:,0]
    return array

def createRanges(u,l,s):
    l1=np.array(range(u,l,s))
    l2=-l1
    ranges=np.array([[l1[i], l2[i]] for i in range(len(l1))])
    r2=np.reshape(ranges, len(l1)*2)
    return r2

def getResolution(f):
    '''Returns x,y,z resolution in mm.'''
    res=list_values(read_line('VisuCoreResolution=',f))
    if len(res)==0:
        try:
            res=list_values(read_line('VisuCoreResolutionMS=',f))
        except AttributeError:
            pass
    return res
    
def checkFileType(file):
    '''Checks the current folder to see what extension the input file has.
    Returns this extension so that it can be easily appended to your filename.

    Inputs:
        file: A string without extension.
    Outputs:
        ext: A string containing the file's extension.'''
    files=glob.glob(file+'.*')
    if file + '.img' in files:
        return '.img'
    elif file + '.nii.gz' in files:
        return '.nii.gz'
    elif file + '.nii' in files:
        return '.nii'

        
