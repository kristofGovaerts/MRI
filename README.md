# MRI
Useful tools for preclinical MRI analysis. Note that all these tools were written in Python 2.7 some time ago - I have no intention of updating them anytime soon. 

## ProcessImages

The file to run here is titled **ProcessImages.py**. It is dependent on several other files in the same folder, such as **ImageTypes.py** and **ImProcessing.py**. The function of this program is to scan a directory containing preclinical MRI scans, analyze the scan type for each, and run an appropriate processing pipeline, such as T2 mapping, perfusion processing, or DTI calculation.

## Conversion

Contains some useful tools for converting data from Bruker to Nifti as well as for sorting scans. 
