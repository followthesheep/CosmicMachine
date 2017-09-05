import numpy as np
import matplotlib as plt
from astropy.io import fits
import os

def mkmask(darklist,indir=None,threshold=0.2,outfile='badpixelmask.fits',medianfile=None):
    '''
    Make a bad pixel mask from a series of darks

    Input
    -----
    darklist - a text file with a list of fits locations

    Keywords
    --------
    threshold - the threshold for a bad pixel

    History
    -------
    2017-07-18 - T. Do 
    
    '''
    infiles = np.loadtxt(darklist,dtype=bytes).astype(str)

    for i in np.arange(len(infiles)):
        if indir is None:
            tmphdu = fits.open(infiles[i])
        else:
            s = os.path.split(infiles[i])
            print(s)
            tmpfile = os.path.join(indir,s[1])
            tmphdu = fits.open(tmpfile)

        data = tmphdu[0].data
        if i == 0:
            stack = np.zeros((data.shape[0],data.shape[1],len(infiles)))
        stack[:,:,i] = data

    # need to set the type in order to match what the pipeline expects
    outarr = np.zeros(data.shape,dtype='uint8')+9  

    medarr = np.median(stack,axis=2)

    badpix = np.where(medarr > threshold)
    print("There are %i bad pixels" % (len(badpix[0])))
    outarr[badpix[0],badpix[1]] = 0
    
    hdu = fits.PrimaryHDU(outarr)
    hdulist = fits.HDUList([hdu])
    print("writing file: "+outfile)
    hdulist.writeto(outfile,overwrite=True)

    if medianfile is not None:
        hdu = fits.PrimaryHDU(medarr)
        
        hdulist = fits.HDUList([hdu,fits.ImageHDU(medarr),fits.ImageHDU(tmphdu[2].data)])
        print("writing file: "+medianfile)
        hdulist.writeto(medianfile,overwrite=True)
        
            
def apply_mask(infile,maskfile,outfile):
    '''
    Apply a bad pixel mask to the third extension of a file

    Inputs
    ------
    infile - input file name of raw
    maskfile - mask file name
    outfile - output file name, can be same as input (will
    overwrite)

    '''

    if os.path.exists(infile):
        hdu = fits.open(infile)
        mask = fits.getdata(maskfile)
        hdu[2].data = mask
        print("writing: "+outfile)
        hdu.writeto(outfile,overwrite=True)
