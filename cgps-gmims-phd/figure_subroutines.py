#import astropy.io.fits
#import sys
import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib import pylab
from astropy.io import fits
#from astropy.wcs import WCS
#import os

def make_axis_lists_3D(hdr):
    
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    nz = hdr['NAXIS3']
   
    dx = hdr['CDELT1']
    dy = hdr['CDELT2']
    dz = hdr['CDELT3']
    
    xpix = hdr['CRPIX1']
    ypix = hdr['CRPIX2']
    zpix = hdr['CRPIX3']
    
    xval = hdr['CRVAL1']
    yval = hdr['CRVAL2']
    zval = hdr['CRVAL3']
    
    x = np.arange(nx)+1-xpix
    y = np.arange(ny)+1-ypix
    z = np.arange(nz)+1-zpix

    lon_ax = x*dx+xval
    lat_ax = y*dy+yval
    z_ax = z*dz+zval
    
    return lon_ax, lat_ax, z_ax

def make_axis_lists_2D(hdr):
    
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
   
    dx = hdr['CDELT1']
    dy = hdr['CDELT2']
    
    xpix = hdr['CRPIX1']
    ypix = hdr['CRPIX2']
    
    xval = hdr['CRVAL1']
    yval = hdr['CRVAL2']
    
    x = np.arange(nx)+1-xpix
    y = np.arange(ny)+1-ypix

    lon_ax = x*dx+xval
    lat_ax = y*dy+yval
    
    return lon_ax, lat_ax
