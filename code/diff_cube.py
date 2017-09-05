import numpy as np
import pylab as plt
from astropy.io import fits
import os

def diff_mask_cube():
    # take a difference between using a bad pixel mask and not using one
    masked_cube = fits.getdata('../reduced/s170517_a018002_Kn3_035_masked.fits')
    nomask_cube = fits.getdata('../reduced/s170517_a018002_Kn3_035_nomask.fits')
    masked_cube = np.transpose(masked_cube,(1,0,2))
    nomask_cube = np.transpose(nomask_cube,(1,0,2))
    print(masked_cube.shape)
    channels = [0,10,100,300]
    plt.clf()
    for i in np.arange(len(channels)):
        plt.subplot(2,2,i+1)
        diff_channel = nomask_cube[:,:,channels[i]]-masked_cube[:,:,channels[i]]
        plt.imshow(diff_channel,origin='lower',vmin=-0.5,vmax=0.5)
        plt.title('Channel %i' % (channels[i]))
        plt.colorbar()
        # look for anomalous points between the two data sets
        bad = np.where(np.abs(diff_channel) > 2.0*np.std(diff_channel))
        print('%i bad points identified' % (len(bad[0])))
        plt.plot(bad[1],bad[0],'o',mfc='None',color='r')
        plt.xlabel('X Pixel')
        plt.ylabel('Y Pixel')
