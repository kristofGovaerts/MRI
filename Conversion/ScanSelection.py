#Script made by Kristof Govaerts(2013), and intended to select only those scans you want for further processing. It leaves original files untouched.

import os
import sys
import glob
import shutil
import Tkinter, tkFileDialog
root = Tkinter.Tk()

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

print "Enter the directory where the conversion results(.hdr, .img, .txt and .mat) are stored."

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
    
def insert(original, new, pos):
    '''Inserts new inside original at pos.'''
    return original[:pos] + new + original[pos:]    

while True:
    rootdir = tkFileDialog.askdirectory(initialdir="/",title='Please select a directory')
    if os.path.isdir(rootdir) is True: #Checks if entered dir exists
        os.chdir(rootdir)
        root.destroy()
        break
    else:
        print "Pathname invalid. Try again."
        continue

txts = glob.glob('*.txt') #list of all txt files

scanlist = list_scans()

numscanlist = []
for scan in scanlist:
    ln = scanlist.index(scan)
    num = str(ln + 1) + '. '
    scan = insert(scan, num, 0)
    numscanlist.append(scan)

jscanlist = '\n'.join(numscanlist) #for readability & index

yn = 'y'

while yn == 'y':
    print "Different scan types found in directory:"
    print jscanlist

    print "Enter the number of the scan type you want to copy. If you want to copy multiple scan types, separate them by spaces (no commas or anything else)."
    numlist = raw_input("> ").split()
    try:
        numlist[:] = [int(x) - 1 for x in numlist]
    except TypeError:
        print "Please enter a number."
        continue
    except NameError:
        print "Please enter a number."
        continue
    except ValueError:
        print "Please enter a number."
        continue

    if max(numlist) + 1 > len(scanlist) or min(numlist) < 0:
        print "Some values in list are invalid."
        continue
    else:
        pass

    for scan in numlist:
        tardir = os.path.join(rootdir, scanlist[scan])
        print "Copying all %s scans to directory %s." % (scanlist[scan], tardir)
        if os.path.isdir(tardir) is False:
            try:
                os.mkdir(tardir)
            except WindowsError:
                tardir=tardir[:-1]
                os.mkdir(tardir)
        else:
            pass

        for txt in txts:
            o = open(txt, 'r')
            c = 'VisuAcquisitionProtocol=' + scanlist[scan] + '\n'
            if c in o: #copies text, img and hdr files to new dir
                shutil.copy2(txt, tardir)
                
                for ext in ('.img','.hdr','.nii'):
                    try:
                        f = txt[:-4] + ext
                        shutil.copy2(f, tardir)
                    except IOError:
                        continue
                o.close()

            else:
                o.close()
    print "All files copied."
    yn = raw_input("Would you like to choose another scan type? y/n\n> ")

    if 'y' in yn:
        pass

    else:
        print "Goodbye!"
        sys.exit()

#Copyright Kristof Govaerts, 2013