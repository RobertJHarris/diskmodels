# test script for reading in fits random group (i.e., vis) files
import sys

try:
    print(sys.version)
    if(sys.version[0] == 2):    
        pyfits_path = '/opt/casa-release-4.6.0-el6/lib/python2.7/site-packages/' # couldn't get the path to work anyway else....
    else:
        pyfits_path = '/home/rjharris//.local/lib/python3.4/site-packages/' # just use mine for now. python 3 v. wasn't installed on eeyore...
    sys.path.append(pyfits_path)
    import pyfits
except:
    print("No pyfits on this system. Check python path / installation.")
    exit(-1)

print("Pyfits found.")


print("astropy found.")
import pprint 
import numpy as np
import sys
import os
try:
    import psutil     
    HAVE_PSUTIL = True
except: 
    print("No psutil on this system. Memory check will not work. [Not a fatal error]")

CC   = 299792458.0
HIGH = 5 
debug = HIGH

    
def make_fits_model(filename,array,pix_scale,refpix):
    os.system('rm ' + filename)
    #print("Writing to blag")
    hdu = pyfits.PrimaryHDU(array)
    hdulist = pyfits.HDUList([hdu])
    prihdr = hdulist[0].header
    print(pix_scale)
    print(refpix)
#    print('the header reads')
    
    #hdulist[0].header
#    prihdr.set('CDELT1', pix_scale)
    # hdulist[0].header
#    prihdr.set('CDELT2', pix_scale)

#prihdr.set('CDELT2', pix_scale)
    #prihdr.set('CRPIX1', refpix)
    #prihdr.set('CRPIX2', refpix)
    hdulist.writeto(filename)
    print("Data written...")
#    data,hdr = pyfits.getdata(filename,1,header=True)
#    print(hdr)
    #hdr['CTYPE1']= 'RA---TAN'
   # hdr['CTYPE2']= 'DEC--TAN'
   # pyfits.writeto(filename,data,hdr)
#   return 
def main():
    file = './data/IRAS16293_Band9.fixed.rebin.rest.uvfits'
    blah = load_data(file)

    return blah

