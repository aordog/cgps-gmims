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

def make_axis_lists_3D_WHAM(hdr):
    
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    nz = hdr['NAXIS3']
   
    dx = hdr['CD1_1']
    dy = hdr['CD2_2']
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


def splice_disk(imagefile,bmin,bmax,lon,lat):
    
    l1 = 0.
    l2 = 180.
    l3 = -179.5
    l4 = -0.5
    
    wb1 = np.where(abs(lat-bmin)<0.25)[0][0]
    wb2 = np.where(abs(lat-bmax)<0.25)[0][0]+1
    
    wl1 = np.where(abs(lon-l1)<0.25)[0][0]+1
    wl2 = np.where(abs(lon-l2)<0.25)[0][0]
    wl3 = np.where(abs(lon-l3)<0.25)[0][0]+1
    wl4 = np.where(abs(lon-l4)<0.25)[0][0]
    
    #print(lat[wb1:wb2])
    #print(lon[wl2:wl1])
    #print(lon[wl4:wl3])
    
    mapleft = imagefile[wb1:wb2,wl4:wl3]
    #print(mapleft.shape)    
    mapright = imagefile[wb1:wb2,wl2:wl1]
    #print(mapright.shape)    
    mapall = np.concatenate((mapleft,mapright),axis=1)
    #print(mapall.shape)
    
    lonleft = lon[wl4:wl3]+360.
    lonright = lon[wl2:wl1]
    lon_disk = np.concatenate((lonleft,lonright))
    #print(lon_disk)
    lat_disk = lat[wb1:wb2]
    
    return(mapall,lon_disk,lat_disk)


def slice_rect_CGPS_v2(imagefile,lminslc,lmaxslc,bminslc,bmaxslc,lon,lat,dG,dC):
    
    numxplt = int((lmaxslc-lminslc)/dG)+1
    numyplt = int((bmaxslc-bminslc)/dG)+1
    
    lon_sub = np.zeros(numxplt)
    lat_sub = np.zeros(numyplt)    
    rectpiece = np.zeros((numyplt,numxplt))
    print(rectpiece.shape)
    
    for j in range(numyplt):
        print(j)
        centb = bminslc+j*dG
        lat_sub[j] = centb
        
        for i in range(numxplt):
            
            centl = lmaxslc-i*dG
            lon_sub[i] = centl
            
            left = centl+dG/2.
            right = centl-dG/2.
            bot = centb-dG/2.
            top = centb+dG/2.
           
            wx = np.where((lon > right) & (lon < left))
            wy = np.where((lat > bot) & (lat < top))
            #print(right,left,top,bot,np.size(wx),np.size(wy))
            tester = imagefile[np.min(wy):np.max(wy),np.min(wx):np.max(wx)].flatten()
            wbad = np.where(np.isnan(tester))
            good = np.size(tester)-np.size(wbad)
            
            if (good > 100):
                rectpiece[j,i] = np.nanmean(imagefile[np.min(wy):np.max(wy),np.min(wx):np.max(wx)])
            else:
                rectpiece[j,i] = np.nan
            #print(rectpiece[j,i])
  
    return(lon_sub,lat_sub,rectpiece)

def zero_moment(hdr,data,peakPI_data,peakthresh,maxFD,perc,lon_arr,lat_arr,RM_arr,dFD):
    
    print('Calculating zeroth moment...')
    
    zero_mom_hdr = hdr.copy()    
    zero_mom_hdr['NAXIS'] = 2    
    del zero_mom_hdr['CRVAL3']
    del zero_mom_hdr['CTYPE3']
    del zero_mom_hdr['CDELT3']
    del zero_mom_hdr['CUNIT3']
    del zero_mom_hdr['NAXIS3']
    del zero_mom_hdr['LAMSQ0']
    
    xsize = data.shape[2]
    ysize = data.shape[1]
    zero_mom_data = np.empty((ysize,xsize))
    print(xsize,ysize)
        
    for i in range(0,xsize):
        print(i)
        #if i==int(xsize/2):
        #    print('halfway')
        for j in range(0,ysize): 
        #for j in range(0,350):
            #print(j)
            wgood = np.where((data[:,j,i] >= peakthresh) & 
                             (data[:,j,i] >= 0.15*peakPI_data[j,i]) &
                             (abs(RM_arr) <= maxFD))
            #print(lon_arr[i],lat_arr[j],np.size(wgood))
            #print(RM_arr[wgood])
            #print('')
            if np.size(wgood) != 0:
                zero_mom_data[j,i] = dFD*np.sum(data[wgood,j,i])
            else:
                zero_mom_data[j,i] = np.nan
        #print('')        
    
    print('done')
    print('')
    return(zero_mom_data,zero_mom_hdr)


def first_moment(hdr,data,peakPI_data,peakthresh,maxFD,perc,lon_arr,lat_arr,RM_arr,dFD,zero_mom):
    
    print('Calculating first moment...')
    
    first_mom_hdr = hdr.copy()    
    first_mom_hdr['NAXIS'] = 2    
    del first_mom_hdr['CRVAL3']
    del first_mom_hdr['CTYPE3']
    del first_mom_hdr['CDELT3']
    del first_mom_hdr['CUNIT3']
    del first_mom_hdr['NAXIS3']
    del first_mom_hdr['LAMSQ0']
    
    xsize = data.shape[2]
    ysize = data.shape[1]
    first_mom_data = np.empty((ysize,xsize))
    print(xsize,ysize)
        
    for i in range(0,xsize):
        #if i==int(xsize/2):
        #    print('halfway')
        for j in range(0,ysize): 
            wgood = np.where((data[:,j,i] >= peakthresh) & 
                             (data[:,j,i] >= 0.15*peakPI_data[j,i]) &
                             (abs(RM_arr) <= maxFD))
            #print(lon_arr[i],lat_arr[j],np.size(wgood))
            #print(RM_arr[wgood])
            #print('')
            if np.size(wgood) != 0:
                first_mom_data[j,i] = dFD*np.sum(data[wgood,j,i]*RM_arr[wgood])/zero_mom[j,i]
            else:
                first_mom_data[j,i] = np.nan
    
    print('done')
    print('')                    
    return(first_mom_data,first_mom_hdr)


def second_moment(hdr,data,peakPI_data,peakthresh,maxFD,perc,lon_arr,lat_arr,RM_arr,dFD,zero_mom,first_mom):
    
    print('Calculating second moment...')
    
    second_mom_hdr = hdr.copy()    
    second_mom_hdr['NAXIS'] = 2    
    del second_mom_hdr['CRVAL3']
    del second_mom_hdr['CTYPE3']
    del second_mom_hdr['CDELT3']
    del second_mom_hdr['CUNIT3']
    del second_mom_hdr['NAXIS3']
    del second_mom_hdr['LAMSQ0']
    
    xsize = data.shape[2]
    ysize = data.shape[1]
    second_mom_data = np.empty((ysize,xsize))
    print(xsize,ysize)
        
    for i in range(0,xsize):
        if i==int(xsize/2):
            print('halfway')
        for j in range(0,ysize): 
            wgood = np.where((data[:,j,i] >= peakthresh) & 
                             (data[:,j,i] >= 0.15*peakPI_data[j,i]) &
                             (abs(RM_arr) <= maxFD))
            #print(lon_arr[i],lat_arr[j],np.size(wgood))
            #print(RM_arr[wgood])
            #print('')
            if np.size(wgood) != 0:
                second_mom_data[j,i] = dFD*np.sum( data[wgood,j,i] * (RM_arr[wgood] - first_mom[j,i])**2 )/zero_mom[j,i]
            else:
                second_mom_data[j,i] = np.nan
    
    print('done')
    print('')                    
    return(second_mom_data,second_mom_hdr)



def Taylor_readin(cat,num):

    EG_lon = []
    EG_lat = []
    EG_RM = []

    for i in range(0,num):
        #print(i)
        if cat[i][46] == '-':
            # negative 3 digit
            lon1 = float(cat[i][47])
            lon2 = float(cat[i][48])
            lon3 = float(cat[i][49])
            lon4 = float(cat[i][51])
            lon5 = float(cat[i][52])
            lon6 = float(cat[i][53])
            lon7 = float(cat[i][54])
            EG_lon.append(-(lon1*100+lon2*10+lon3+lon4*0.1+lon5*0.01+lon6*0.001+lon7*0.0001))
        else:
            if cat[i][47] != ' ':
                if cat[i][47] == '-':
                    # negative 2 digit            
                    lon2 = float(cat[i][48])
                    lon3 = float(cat[i][49])
                    lon4 = float(cat[i][51])
                    lon5 = float(cat[i][52])
                    lon6 = float(cat[i][53])
                    lon7 = float(cat[i][54])
                    EG_lon.append(-(lon2*10+lon3+lon4*0.1+lon5*0.01+lon6*0.001+lon7*0.0001))
                else:
                    # positive 3 digit
                    lon1 = float(cat[i][47])
                    lon2 = float(cat[i][48])
                    lon3 = float(cat[i][49])
                    lon4 = float(cat[i][51])
                    lon5 = float(cat[i][52])
                    lon6 = float(cat[i][53])
                    lon7 = float(cat[i][54])
                    EG_lon.append(lon1*100+lon2*10+lon3+lon4*0.1+lon5*0.01+lon6*0.001+lon7*0.0001)
            else:
                if cat[i][48] != ' ':
                    if cat[i][48] == '-':
                        # negative 1 digit            
                        lon3 = float(cat[i][49])
                        lon4 = float(cat[i][51])
                        lon5 = float(cat[i][52])
                        lon6 = float(cat[i][53])
                        lon7 = float(cat[i][54])
                        EG_lon.append(-(lon3+lon4*0.1+lon5*0.01+lon6*0.001+lon7*0.0001))
                    else:
                        # positive 2 digit
                        lon2 = float(cat[i][48])
                        lon3 = float(cat[i][49])
                        lon4 = float(cat[i][51])
                        lon5 = float(cat[i][52])
                        lon6 = float(cat[i][53])
                        lon7 = float(cat[i][54])
                        EG_lon.append(lon2*10+lon3+lon4*0.1+lon5*0.01+lon6*0.001+lon7*0.0001)
                else:
                    if cat[i][49] != ' ':
                        # positive 1 digit
                        lon3 = float(cat[i][49])
                        lon4 = float(cat[i][51])
                        lon5 = float(cat[i][52])
                        lon6 = float(cat[i][53])
                        lon7 = float(cat[i][54])
                        EG_lon.append(lon3+lon4*0.1+lon5*0.01+lon6*0.001+lon7*0.0001)
    

    
    for i in range(0,num):   
        if cat[i][58] == '-':
            # negative 2 digit
            lat1 = float(cat[i][59])
            lat2 = float(cat[i][60])
            lat3 = float(cat[i][62])
            lat4 = float(cat[i][63])
            lat5 = float(cat[i][64])
            lat6 = float(cat[i][65])
            EG_lat.append(-(lat1*10+lat2+lat3*0.1+lat4*0.01+lat5*0.001+lat6*0.0001))
        else:
            if cat[i][59] != ' ':
                if cat[i][59] == '-':
                    # negative 1 digit            
                    lat2 = float(cat[i][60])
                    lat3 = float(cat[i][62])
                    lat4 = float(cat[i][63])
                    lat5 = float(cat[i][64])
                    lat6 = float(cat[i][65])
                    EG_lat.append(-(lat2+lat3*0.1+lat4*0.01+lat5*0.001+lat6*0.0001))
                else:
                    # positive 2 digit
                    lat1 = float(cat[i][59])
                    lat2 = float(cat[i][60])
                    lat3 = float(cat[i][62])
                    lat4 = float(cat[i][63])
                    lat5 = float(cat[i][64])
                    lat6 = float(cat[i][65])
                    EG_lat.append(lat1*10+lat2+lat3*0.1+lat4*0.01+lat5*0.001+lat6*0.0001)
            else:
                if cat[i][60] != ' ':
                    # positive 1 digit
                    lat2 = float(cat[i][60])
                    lat3 = float(cat[i][62])
                    lat4 = float(cat[i][63])
                    lat5 = float(cat[i][64])
                    lat6 = float(cat[i][65])
                    EG_lat.append(lat2+lat3*0.1+lat4*0.01+lat5*0.001+lat6*0.0001)

                    
    for i in range(0,num):   
        if cat[i][126] == '-':
            # negative 4 digit
            RM1 = float(cat[i][127])
            RM2 = float(cat[i][128])
            RM3 = float(cat[i][129])
            RM4 = float(cat[i][130])
            RM5 = float(cat[i][132])
            EG_RM.append(-(RM1*1000+RM2*100+RM3*10+RM4+RM5*0.1))
            #EG_lon.append(-(lon1*100+lon2*10+lon3+lon4*0.1+lon5*0.01+lon6*0.001+lon7*0.0001))
        else:
            if cat[i][127] != ' ':
                if cat[i][127] == '-':
                    # negative 3 digit            
                    RM2 = float(cat[i][128])
                    RM3 = float(cat[i][129])
                    RM4 = float(cat[i][130])
                    RM5 = float(cat[i][132])
                    EG_RM.append(-(RM2*100+RM3*10+RM4+RM5*0.1))
                else:
                    # positive 4 digit
                    RM1 = float(cat[i][127])
                    RM2 = float(cat[i][128])
                    RM3 = float(cat[i][129])
                    RM4 = float(cat[i][130])
                    RM5 = float(cat[i][132])
                    EG_RM.append(RM1*1000+RM2*100+RM3*10+RM4+RM5*0.1)
            else:
                if cat[i][128] != ' ':
                    if cat[i][128] == '-':
                        # negative 2 digit            
                        RM3 = float(cat[i][129])
                        RM4 = float(cat[i][130])
                        RM5 = float(cat[i][132])
                        EG_RM.append(-(RM3*10+RM4+RM5*0.1))
                    else:
                        # positive 3 digit
                        RM2 = float(cat[i][128])
                        RM3 = float(cat[i][129])
                        RM4 = float(cat[i][130])
                        RM5 = float(cat[i][132])
                        EG_RM.append(RM2*100+RM3*10+RM4+RM5*0.1)
                else:
                    if cat[i][129] != ' ':
                        if cat[i][129] == '-':
                            # negative 1 digit            
                            RM4 = float(cat[i][130])
                            RM5 = float(cat[i][132])
                            EG_RM.append(-(RM4+RM5*0.1))
                        else:
                            # positive 2 digit
                            RM3 = float(cat[i][129])
                            RM4 = float(cat[i][130])
                            RM5 = float(cat[i][132])
                            EG_RM.append(RM3*10+RM4+RM5*0.1)
                    else:
                        if cat[i][130] != ' ':
                            # positive 1 digit
                            RM4 = float(cat[i][130])
                            RM5 = float(cat[i][132])
                            EG_RM.append(RM4+RM5*0.1)

    EG_lon_out = np.asarray(EG_lon) 
    EG_lat_out = np.asarray(EG_lat)
    EG_RM_out = np.asarray(EG_RM)
                                       
    return(EG_lon_out,EG_lat_out,EG_RM_out)




def zero_moment_v2(hdr,data,peakPI_data,peakthresh,maxFD,perc,lon_arr,lat_arr,RM_arr,dFD):
    
    print('Calculating zeroth moment...')
    
    zero_mom_hdr = hdr.copy()    
    zero_mom_hdr['NAXIS'] = 2    
    del zero_mom_hdr['CRVAL3']
    del zero_mom_hdr['CTYPE3']
    del zero_mom_hdr['CDELT3']
    del zero_mom_hdr['CUNIT3']
    del zero_mom_hdr['NAXIS3']
    del zero_mom_hdr['LAMSQ0']
    
    xsize = data.shape[2]
    ysize = data.shape[1]
    zero_mom_data = np.empty((ysize,xsize))
    zero_mom_data[:,:] = np.nan
    print(xsize,ysize)
    
    data_use = data.copy()
    wbigFD = np.where(abs(RM_arr) > maxFD)
    data_use[wbigFD,:,:] = np.nan
    
    for i in range(0,xsize):
        if i==int(xsize/2):
            print('halfway')
        for j in range(0,ysize): 
            
            comparer = max(peakthresh,0.15*peakPI_data[j,i])
            dataSlice = data_use[:,j,i]
            zero_mom_data[j,i] = np.nansum(dFD*dataSlice[dataSlice >= comparer])

    
    print('done')
    print('')
    return(zero_mom_data,zero_mom_hdr)


def first_moment_v2(hdr,data,peakPI_data,peakthresh,maxFD,perc,lon_arr,lat_arr,RM_arr,dFD,zero_mom):
    
    print('Calculating first moment...')
    
    first_mom_hdr = hdr.copy()    
    first_mom_hdr['NAXIS'] = 2    
    del first_mom_hdr['CRVAL3']
    del first_mom_hdr['CTYPE3']
    del first_mom_hdr['CDELT3']
    del first_mom_hdr['CUNIT3']
    del first_mom_hdr['NAXIS3']
    del first_mom_hdr['LAMSQ0']
    
    xsize = data.shape[2]
    ysize = data.shape[1]
    first_mom_data = np.empty((ysize,xsize))
    first_mom_data[:,:] = np.nan
    print(xsize,ysize)
    
    data_use = data.copy()
    wbigFD = np.where(abs(RM_arr) > maxFD)
    data_use[wbigFD,:,:] = np.nan
        
    for i in range(0,xsize):
        if i==int(xsize/2):
            print('halfway')
        for j in range(0,ysize): 
            comparer = max(peakthresh,0.15*peakPI_data[j,i])
            dataSlice = data_use[:,j,i]
            dataSliceRM = dataSlice*RM_arr
            first_mom_data[j,i] = dFD*dataSliceRM[dataSlice >= comparer].sum()/zero_mom[j,i]

    
    print('done')
    print('')                    
    return(first_mom_data,first_mom_hdr)


