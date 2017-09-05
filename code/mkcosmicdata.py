import numpy as np
import pylab as plt
from astropy.io import fits
from astropy.stats import sigma_clip
import os


def mk_median_sky(outmedian='../data/skies/skies_median.fits',mask='../mask/badpixelmask2017.fits.gz'):
    # identify cosmic rays in a stack of skies
    directory = '../data/skies/'
    files = ['s170517_a018002.fits','s170518_a016003.fits','s170519_a011002.fits',
             's170719_a009002.fits','s170727_a010002.fits']


    for i in np.arange(len(files)):
        tmpfile = os.path.join(directory,files[i])
        hdu = fits.open(tmpfile)
        data = hdu[0].data
        if i == 0:
            s = data.shape
            stack = np.zeros((s[0],s[1],len(files)))
        stack[:,:,i] = data

    med = np.median(stack, axis=2)

    bpm = fits.getdata(mask)

    medhdu = fits.PrimaryHDU(med)
    medhdu.header = hdu[0].header # use the header from the last image for the median
    hdulist = fits.HDUList([medhdu,fits.ImageHDU(hdu[1].data),fits.ImageHDU(bpm)])
    hdulist.writeto(outmedian)

def id_cosmic(mediansky='../data/skies/skies_median.fits',mask='../mask/badpixelmask2017.fits.gz'):
    # identify cosmic rays in a stack of skies
    directory = '../data/skies/'
    files = ['s170517_a018002.fits','s170518_a016003.fits','s170519_a011002.fits',
             's170719_a009002.fits','s170727_a010002.fits']

    median_sky = fits.getdata(mediansky)
    print("identifying cosmic rays")
    for i in np.arange(len(files)):
        tmpfile = os.path.join(directory,files[i])
        hdu = fits.open(tmpfile)
        data = hdu[0].data

        diff = data - median_sky
        outhdu = fits.PrimaryHDU(diff)
        outhdu.data = diff
        parts = os.path.splitext(files[i])
        outhdu.writeto(parts[0]+'_diff.fits',overwrite=True)

        s=sigma_clip(diff,sigma=4,sigma_lower=50,sigma_upper=5)
        good_arr = np.zeros(np.shape(diff))
        good_arr[s.mask] = 7 # give cosmic rays a value of 7

        cosmichdu = fits.PrimaryHDU(good_arr)
        cosmichdu.writeto(parts[0]+'_cosmic.fits',overwrite=True)

def create_cosmic_cube(infile,cosmicfile=None,npix=5,outfile=None):
    '''
    create a data cube centered around each identified cosmic ray point

    Inputs
    ------
    infile - inputfile of the data to extract cosmic rays from

    Keywords
    --------
    cosmicfile - the fits file that contains the cosmic ray mask (pixels labeled as 7)
    npix - the number of pixels to extract around each cosmic ray
    outfile - output file name of the stacks (default: input file +  _cosmic_cube.fits)
    '''
    hdu = fits.open(infile)
    data = hdu[0].data

    # assume the mask is in the second extension if not specifically called
    if cosmicfile is None:
        cosmic_mask = hdu[2].data
    else:
        cosmic_mask = fits.getdata(cosmicfile)

    # cosmic rays are labeled with a 7
    cosmic_ind = np.where(cosmic_mask == 7)
    outarr = np.zeros((npix,npix,len(cosmic_ind[0])))

    if (npix % 2) == 0:
        extra = 0
    else:
        extra = 1

    data_shape = np.shape(data)
    cosmic_bool = np.ones((len(cosmic_ind[0])),dtype=bool) # this array represents whether the slice is a valid cosmic ray

    for i in np.arange(len(cosmic_ind[0])):
        pos = [cosmic_ind[0][i],cosmic_ind[1][i]]
        if (pos[0] > int(npix/2)) & (pos[0] < data_shape[0]-int(npix/2)+extra) & \
           (pos[1] > int(npix/2)) & (pos[1] < data_shape[1]-int(npix/2)+extra):
                subarr = data[pos[0]-int(npix/2):pos[0]+int(npix/2)+extra,
                              pos[1]-int(npix/2):pos[1]+int(npix/2)+extra]
                outarr[:,:,i] = subarr
                cosmic_bool[i] = True
        else:
            cosmic_bool[i] = False

    # reindex the array to remove the slices without cosmic rays
    outarr = outarr[:,:,cosmic_bool]

    outhdu = fits.PrimaryHDU(outarr)
    outhdu.data = outarr
    parts = os.path.splitext(infile)
    if outfile is None:
        outhdu.writeto(parts[0]+'_cosmic_cube.fits',overwrite=True)
    else:
        outhdu.writeto(outfile,overwrite=True)

def prep_input(data,positions,npix=3):
    '''
    Prepare a dataset for processing by flattening the arrays and returning
    an array of the form: [nsamples, nfeatures]

    Currently will not work near edges

    Inputs
    ------
    data - input 2D array of data
    positions - central pixel of the data [npoints, positions]

    Keywords
    ------
    npix - the box size for extraction (default: 3)

    '''
    s = np.shape(positions)
    data_shape = np.shape(data)
    nfeatures = npix*npix
    outarr = np.ones((s[0],nfeatures))

    # need to account for odd vs even npix for subscripting
    if (npix % 2) == 0:
        extra = 0
    else:
        extra = 1

    for i in np.arange(s[0]):
        pos = positions[i,:]
        pos = np.flip(pos,axis=0) # this takes care of the weird python flip of arrays
        if (pos[0] > int(npix/2)) & (pos[0] < data_shape[0]-int(npix/2)+extra) & \
           (pos[1] > int(npix/2)) & (pos[1] < data_shape[1]-int(npix/2)+extra):
                subarr = data[pos[0]-int(npix/2):pos[0]+int(npix/2)+extra,
                              pos[1]-int(npix/2):pos[1]+int(npix/2)+extra]

                #subarr = subarr/np.median(subarr)
                outarr[i,:] = subarr.flatten()
    return outarr

def test():
    print('testing 2')
