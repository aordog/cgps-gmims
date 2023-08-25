import os
import subprocess
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import units as u
from ipywidgets import interact
import matplotlib.patches as patches
from importlib import reload 
from ipywidgets import interact, interactive, fixed, interact_manual, Layout
import ipywidgets as widgets
from mpl_point_clicker import clicker
from tqdm import tqdm
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord

def do_RM_compare():
    
    hdu_RM_CG = fits.open('/srv/aordog/cgps_gmims_data/RM_CG_conv4_regrd.fits')
    hdu_RM_G  = fits.open('/srv/aordog/cgps_gmims_data/RM_G_conv4_regrd.fits')
    RM_CG     = hdu_RM_CG[0].data
    RM_G      = hdu_RM_G[0].data
    RM_CG_hdr = hdu_RM_CG[0].header

    RM_CG[RM_CG==0] = np.nan
    RM_G[RM_G==0] = np.nan

    rvalue_CG = hdu_RM_CG[2].data
    rvalue_G  = hdu_RM_G[2].data
    
    stderr_CG = hdu_RM_CG[4].data
    stderr_G  = hdu_RM_G[4].data

    PAint_CG  = hdu_RM_CG[1].data
    PAint_G   = hdu_RM_G[1].data

    hdu_PI_CG = fits.open('/srv/aordog/cgps_gmims_data/PI_CG_conv4_regrd_avg_PI.fits')
    hdu_PI_G = fits.open('/srv/aordog/cgps_gmims_data/PI_G_conv4_regrd_avg_PI.fits')
    PI_CG     = hdu_PI_CG[0].data
    PI_G     = hdu_PI_G[0].data

    hdu_PA_A_CG = fits.open('/srv/aordog/cgps_gmims_data/PA_A_CG_conv4_regrd.fits')
    hdu_PA_B_CG = fits.open('/srv/aordog/cgps_gmims_data/PA_B_CG_conv4_regrd.fits')
    hdu_PA_C_CG = fits.open('/srv/aordog/cgps_gmims_data/PA_C_CG_conv4_regrd.fits')
    hdu_PA_D_CG = fits.open('/srv/aordog/cgps_gmims_data/PA_D_CG_conv4_regrd.fits')

    hdu_PA_A_G = fits.open('/srv/aordog/cgps_gmims_data/PA_A_G_conv4_regrd.fits')
    hdu_PA_B_G = fits.open('/srv/aordog/cgps_gmims_data/PA_B_G_conv4_regrd.fits')
    hdu_PA_C_G = fits.open('/srv/aordog/cgps_gmims_data/PA_C_G_conv4_regrd.fits')
    hdu_PA_D_G = fits.open('/srv/aordog/cgps_gmims_data/PA_D_G_conv4_regrd.fits')

    PA_A_CG = hdu_PA_A_CG[0].data
    PA_B_CG = hdu_PA_B_CG[0].data
    PA_C_CG = hdu_PA_C_CG[0].data
    PA_D_CG = hdu_PA_D_CG[0].data

    PA_A_G = hdu_PA_A_G[0].data
    PA_B_G = hdu_PA_B_G[0].data
    PA_C_G = hdu_PA_C_G[0].data
    PA_D_G = hdu_PA_D_G[0].data
  
    PA_CG_arr = [PA_A_CG,PA_B_CG,PA_C_CG,PA_D_CG]

    make_the_plot(RM_CG, PI_CG, RM_CG_hdr, PA_CG_arr, rvalue_CG, PAint_CG, stderr_CG,
                  datamax1=300,datamax2=1)

    return


def mouse_event(event,ax,PA_arr,rvalue,PAint,data1,hdr,stderr):
        
    print('x: {} and y: {}'.format(event.xdata, event.ydata))
    ii = int(np.round(event.xdata))
    jj = int(np.round(event.ydata))
    print(ii,jj)
    ax.cla()

    freq = np.array([1406.9,1413.8,1427.4,1434.3])
    lbd2 = ((3e8)/(freq*1e6))**2
    lbd2_ext = np.linspace(0.0435,0.0456,100)
    
    PAint_pt = PAint[jj,ii]
    RM_pt = data1[jj,ii]
    PA_arr_plt = np.array([PA_arr[0][jj,ii],PA_arr[1][jj,ii],
                           PA_arr[2][jj,ii],PA_arr[3][jj,ii]])*180/np.pi
    
    ax.plot(lbd2_ext,PAint_pt*180/np.pi+RM_pt*lbd2_ext*180/np.pi,color='blue')
    ax.scatter(lbd2,PA_arr_plt,color='C0')
    ax.scatter(lbd2,PA_arr_plt+180,color='C1')
    ax.scatter(lbd2,PA_arr_plt-180,color='C2')
    ax.set_xlim(0.0435,0.0456)
    ax.set_xticks([0.0435,0.044,0.0445,0.045,0.0455])
    ax.set_ylim(-180,180)
    ax.set_yticks([-180,-135,-90,-45,0,45,90,135,180])
    ax.set_xlabel(r'$\lambda^2$ (m$^2$)')
    ax.set_ylabel('PA (deg.)')
    ax.grid()
    
    point_coords = pixel_to_skycoord(ii,jj,WCS(hdr))
    print(point_coords)
    ax.text(0.0437,-260,'l='+str(np.round(point_coords.l.deg,2))+
                       ' b='+str(np.round(point_coords.b.deg,2)),fontsize=16)
    ax.text(0.0437,-290,'r='+str(np.round(rvalue[jj,ii],2)),fontsize=16)
    ax.text(0.0437,-320,'RM='+str(np.round(RM_pt,2))+r' rad m $^{-2}$',fontsize=16)
    ax.text(0.0437,-350,r'$\sigma_{RM}$='+str(np.round(stderr[jj,ii],1))+r' rad m $^{-2}$',fontsize=16)

    plt.show()

    return 
                                 
def make_the_plot(data1,data2,hdr,PA_arr,rvalue,PAint,stderr,
                  datamax1=100, datamax2=0.5,llim=[82,52], blim=[-7,10],*args,**kwargs):
  
    c = SkyCoord(llim, blim, frame=Galactic, unit="deg")
        
    fig1 = plt.figure(figsize=(20,10))
    
    ax1  = fig1.add_subplot(221, projection=WCS(hdr).celestial)
    im1  = ax1.imshow(data1, origin='lower', vmin=-datamax1, vmax=datamax1,cmap='Spectral_r')
    ax1.set_xlim(WCS(hdr).world_to_pixel(c)[0])
    ax1.set_ylim(WCS(hdr).world_to_pixel(c)[1])
    cbar1 = fig1.colorbar(im1,orientation='horizontal',fraction=0.06)
    #cbar1.set_label(r'RM (rad m $^{-2}$')
    
    ax2  = fig1.add_subplot(223, projection=WCS(hdr).celestial, sharex=ax1, sharey=ax1)
    im2  = ax2.imshow(data2, origin='lower', vmin=0, vmax=datamax2,cmap='cubehelix')
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())
    cbar2 = fig1.colorbar(im2,orientation='horizontal',fraction=0.06)
    cbar2.set_label('PI (K)')
    
    ax  = fig1.add_subplot(222)
        
    cid = fig1.canvas.mpl_connect('button_press_event', 
                                  lambda event: mouse_event(event,ax,PA_arr,rvalue,PAint,data1,hdr,stderr))
    klicker = clicker(ax1, ["s"], markers=["x"])
    ax1.get_legend().remove()
    
    plt.tight_layout()
    plt.show()
            
    return


if __name__ =='__main__': 
    do_RM_compare()