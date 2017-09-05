import numpy as np
import pylab as plt
from astropy.io import fits
from astropy.stats import sigma_clip
import os

# deals with creating valid data for cosmic ray training purposes

def mkvalidcube(infile,cosmicfile=None,mask=None,rectfile='../data/s170508_c008___infl_Kn3_035.fits.gz',
                tracefile='../data/s170508_c008___infl_Kn3_035_trace.npy',narrowband=True,
                nsim=500,npix=5,debug=False,outfile=None):
    '''
    Make a cube of valid data points for training

    Input
    -----
    infile - input filename of raw image

    Keywords
    --------
    mask - bad pixel mask (default: None)
    narrowband - the filter is a narrowband, so should split the data up into
    3 regions
    '''
    hdu = fits.open(infile)
    data = hdu[0].data

    # assume the mask is in the second extension if not specifically called
    if cosmicfile is None:
        cosmic_mask = hdu[2].data
    else:
        cosmic_mask = fits.getdata(cosmicfile)


    if mask is not None:
        bpm = fits.getdata(mask)
    else:
        bpm = np.zeros(data.shape)+9

    totalmask = bpm + cosmic_mask

    # load the rectmat and the trace
    hdu = fits.open(rectfile)
    # the slices for the rect matrix are in the second extension
    inds = hdu[0].data
    matrix = hdu[2].data # should be shaped (1216, 16, 2048)
    s = np.shape(matrix)

    a = np.load(tracefile)
    b = a.item(0)
    data_range = np.zeros(s[2],dtype=bool)
    if narrowband:
        # Select a range in x that works for the entire detector
        data_range[250:550] = True
        data_range[710:1000] = True
        data_range[1200:1450] = True
    else:
        data_range[100:2000] = True

    if debug:
        plt.clf()
        nsim = 9

    # pick a random x position as well as a random slice
    random_x = np.array(np.fix(np.random.uniform(low=npix/2,high=s[2]-npix,size=nsim)),dtype=int)
    random_slice = np.array(np.fix(np.random.uniform(low=0,high=len(b),size=nsim)),dtype=int)
    valid_points = np.ones(nsim,dtype=bool)

    # need to account for odd vs even npix for subscripting
    if (npix % 2) == 0:
        extra = 0
    else:
        extra = 1

    outarr = np.zeros((npix,npix,nsim))


    for i in np.arange(len(random_x)):
        tempslice = 'slice'+str(random_slice[i])
        # find the y position of the point in the random slice of data
        ypos = (b[tempslice][2])[random_x[i]]
        if (ypos > npix/2) & (ypos > 0):
            ypos = int(ypos + inds[random_slice[i],0])

            # check the mask first to see if this is valid
            if totalmask[ypos,random_x[i]] == 9:
                subarr = data[ypos-int(npix/2):ypos+int(npix/2)+extra,
                         random_x[i]-int(npix/2):random_x[i]+int(npix/2)+extra]
                outarr[:,:,i] = subarr
                if debug:
                    plt.subplot(3,3,i+1)
                    plt.imshow(subarr,origin='bottom')

            else:
                valid_points[i] = False
        else:
            valid_points[i] = False

    # reindex the array to remove the slices without cosmic rays

    outarr = outarr[:,:,valid_points]

    outhdu = fits.PrimaryHDU(outarr)
    outhdu.data = outarr
    parts = os.path.splitext(infile)
    if outfile is None:
        outhdu.writeto(parts[0]+'_valid_cube.fits',overwrite=True)
    else:
        outhdu.writeto(outfile,overwrite=True)

def mkohlinecube(infile,cosmicfile=None,mask=None,
                npix=5,debug=False,outfile=None,
                ohfile='../data/OH_positions.dat'):
    '''
    Make a cube of the OH lines locations. Currently, the OH line list is only
    computed for Kbb 2017
    '''
    hdu = fits.open(infile)
    data = hdu[0].data
    data_shape = data.shape
    # assume the mask is in the second extension if not specifically called
    if cosmicfile is None:
        cosmic_mask = hdu[2].data
    else:
        cosmic_mask = fits.getdata(cosmicfile)


    if mask is not None:
        bpm = fits.getdata(mask)
    else:
        bpm = np.zeros(data.shape)

    totalmask = bpm + cosmic_mask

    # need to account for odd vs even npix for subscripting
    if (npix % 2) == 0:
        extra = 0
    else:
        extra = 1

    wave,x,y,flux = np.loadtxt('../data/OH_positions.dat',skiprows=1,unpack=True)
    x = np.array(x,dtype=int)
    y = np.array(y,dtype=int)
    n = len(x)
    print("number of points: "+str(n))
    outarr = np.zeros((npix,npix,n))
    valid_points = np.ones(n,dtype=bool)
    for i in np.arange(n):
        if (x[i] > int(npix/2)) & (x[i] < data_shape[0]-int(npix/2)+extra) & \
           (y[i] > int(npix/2)) & (y[i] < data_shape[1]-int(npix/2)+extra) & \
           (totalmask[y[i],x[i]] > 0):
            subarr = data[y[i]-int(npix/2):y[i]+int(npix/2)+extra,x[i]-int(npix/2):x[i]+int(npix/2)+extra]
            outarr[:,:,i] = subarr
        else:
            valid_points[i] = False
    # reindex the array to remove the slices without cosmic rays
    #print("number of valid points: "+str(len(np.where(valid_points))))
    outarr = outarr[:,:,valid_points]

    outhdu = fits.PrimaryHDU(outarr)
    outhdu.data = outarr
    parts = os.path.splitext(infile)
    if outfile is None:
        outhdu.writeto(parts[0]+'_ohline_cube.fits',overwrite=True)
    else:
        outhdu.writeto(outfile,overwrite=True)
