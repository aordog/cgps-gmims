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

def do_RM_compare(dir_in,fs=12,llim=[82,52], blim=[-7,10],*args,**kwargs):
    
    hdu_RM_CG = fits.open(dir_in+'RM_CG_conv4_regrd.fits')
    hdu_RM_G  = fits.open(dir_in+'RM_G_conv4_regrd.fits')
    hdu_FD_G  = fits.open(dir_in+'phi_peak_regrd.fits')
    #hdu_RM_G  = fits.open(dir_in+'RM_G_regrd.fits')
    hdu_RM_C  = fits.open(dir_in+'RM_C_conv4_regrd.fits')
    RM_CG     = hdu_RM_CG[0].data
    RM_G      = hdu_RM_G[0].data
    FD_G      = hdu_FD_G[0].data
    RM_C      = hdu_RM_C[0].data
    hdr = hdu_RM_CG[0].header

    RM_CG[RM_CG==0] = np.nan
    RM_G[RM_G==0] = np.nan
    FD_G[FD_G==0] = np.nan
    RM_C[RM_C==0] = np.nan

    rvalue_CG = hdu_RM_CG[2].data
    rvalue_G  = hdu_RM_G[2].data
    rvalue_C  = hdu_RM_C[2].data
    
    stderr_CG = hdu_RM_CG[4].data
    stderr_G  = hdu_RM_G[4].data
    stderr_C  = hdu_RM_C[4].data

    PAint_CG  = hdu_RM_CG[1].data
    PAint_G   = hdu_RM_G[1].data
    PAint_C   = hdu_RM_C[1].data

    hdu_PI_CG = fits.open(dir_in+'PI_CG_conv4_regrd_avg_PI.fits')
    #hdu_PI_G  = fits.open(dir_in+'PI_G_conv4_regrd_avg_PI.fits')
    hdu_PI_G  = fits.open(dir_in+'PI_G_conv4_regrd_avg_PI.fits')
    hdu_PI_C  = fits.open(dir_in+'PI_C_conv4_regrd_avg_PI.fits')
    PI_CG     = hdu_PI_CG[0].data
    PI_G      = hdu_PI_G[0].data
    PI_C      = hdu_PI_C[0].data

    hdu_PA_A_CG = fits.open(dir_in+'PA_A_CG_conv4_regrd.fits')
    hdu_PA_B_CG = fits.open(dir_in+'PA_B_CG_conv4_regrd.fits')
    hdu_PA_C_CG = fits.open(dir_in+'PA_C_CG_conv4_regrd.fits')
    hdu_PA_D_CG = fits.open(dir_in+'PA_D_CG_conv4_regrd.fits')

    hdu_PA_A_G = fits.open(dir_in+'PA_A_G_conv4_regrd.fits')
    hdu_PA_B_G = fits.open(dir_in+'PA_B_G_conv4_regrd.fits')
    hdu_PA_C_G = fits.open(dir_in+'PA_C_G_conv4_regrd.fits')
    hdu_PA_D_G = fits.open(dir_in+'PA_D_G_conv4_regrd.fits')
    #hdu_PA_A_G = fits.open(dir_in+'PA_A_G_regrd.fits')
    #hdu_PA_B_G = fits.open(dir_in+'PA_B_G_regrd.fits')
    #hdu_PA_C_G = fits.open(dir_in+'PA_C_G_regrd.fits')
    #hdu_PA_D_G = fits.open(dir_in+'PA_D_G_regrd.fits')
    
    hdu_PA_A_C = fits.open(dir_in+'PA_A_C_conv4_regrd.fits')
    hdu_PA_B_C = fits.open(dir_in+'PA_B_C_conv4_regrd.fits')
    hdu_PA_C_C = fits.open(dir_in+'PA_C_C_conv4_regrd.fits')
    hdu_PA_D_C = fits.open(dir_in+'PA_D_C_conv4_regrd.fits')

    PA_A_CG = hdu_PA_A_CG[0].data
    PA_B_CG = hdu_PA_B_CG[0].data
    PA_C_CG = hdu_PA_C_CG[0].data
    PA_D_CG = hdu_PA_D_CG[0].data

    PA_A_G = hdu_PA_A_G[0].data
    PA_B_G = hdu_PA_B_G[0].data
    PA_C_G = hdu_PA_C_G[0].data
    PA_D_G = hdu_PA_D_G[0].data
    
    PA_A_C = hdu_PA_A_C[0].data
    PA_B_C = hdu_PA_B_C[0].data
    PA_C_C = hdu_PA_C_C[0].data
    PA_D_C = hdu_PA_D_C[0].data
  
    PA_G  = [PA_A_G,PA_B_G,PA_C_G,PA_D_G]
    PA_C  = [PA_A_C,PA_B_C,PA_C_C,PA_D_C]
    PA_CG = [PA_A_CG,PA_B_CG,PA_C_CG,PA_D_CG]
    
    data_list = [FD_G,RM_G,PI_G, RM_C,PI_C, RM_CG,PI_CG]
    r_list = [rvalue_G,rvalue_C,rvalue_CG]
    err_list = [stderr_G,stderr_C,stderr_CG]
    PAint_list = [PAint_G,PAint_C,PAint_CG]

    make_the_plot(data_list,r_list,err_list,PAint_list,hdr,PA_G,PA_C,PA_CG,fs,
                  llim,blim,datamax1=300,datamax2=1)

    return


def mouse_event(event,ax,data_list,r_list,err_list,PAint_list,hdr,axs,
                PA_G,PA_C,PA_CG,fs):
        
    print('x: {} and y: {}'.format(event.xdata, event.ydata))
    ii = int(np.round(event.xdata))
    jj = int(np.round(event.ydata))
    print(ii,jj)
    ax.cla()

    freq = np.array([1406.9,1413.8,1427.4,1434.3])
    lbd2 = ((3e8)/(freq*1e6))**2
    lbd2_ext = np.linspace(0.0435,0.0456,100)
    
    #PAint_pt = PAint[jj,ii]
    RM_pt_G = data_list[1][jj,ii]
    RM_pt_C = data_list[3][jj,ii]
    RM_pt_CG = data_list[5][jj,ii]
    
    PAint_pt_G = PAint_list[0][jj,ii]
    PAint_pt_C = PAint_list[1][jj,ii]
    PAint_pt_CG = PAint_list[2][jj,ii]
    
    PA_arr_C = np.array([PA_C[0][jj,ii],PA_C[1][jj,ii],
                         PA_C[2][jj,ii],PA_C[3][jj,ii]])*180/np.pi
    PA_arr_G = np.array([PA_G[0][jj,ii],PA_G[1][jj,ii],
                         PA_G[2][jj,ii],PA_G[3][jj,ii]])*180/np.pi
    PA_arr_CG= np.array([PA_CG[0][jj,ii],PA_CG[1][jj,ii],
                         PA_CG[2][jj,ii],PA_CG[3][jj,ii]])*180/np.pi
    
    ax.plot(lbd2_ext,PAint_pt_G*180/np.pi+RM_pt_G*lbd2_ext*180/np.pi,color='C0')
    ax.plot(lbd2_ext,PAint_pt_C*180/np.pi+RM_pt_C*lbd2_ext*180/np.pi,color='C1')
    ax.plot(lbd2_ext,PAint_pt_CG*180/np.pi+RM_pt_CG*lbd2_ext*180/np.pi,color='C2')
    ax.scatter(lbd2,PA_arr_G,color='C0',label='GMIMS')
    ax.scatter(lbd2,PA_arr_C,color='C1',label='CGPS')
    ax.scatter(lbd2,PA_arr_CG,color='C2',label='CGPS+GMIMS')
    #ax.scatter(lbd2,PA_arr_plt+180,color='C1')
    #ax.scatter(lbd2,PA_arr_plt-180,color='C2')
    ax.set_xlim(0.0435,0.0456)
    ax.set_xticks([0.0435,0.044,0.0445,0.045,0.0455])
    ax.set_ylim(-180,180)
    ax.set_yticks([-180,-135,-90,-45,0,45,90,135,180])
    ax.set_xlabel(r'$\lambda^2$ (m$^2$)')
    ax.set_ylabel('PA (deg.)')
    ax.grid()
    ax.legend(fontsize=fs-4)
    
    point_coords = pixel_to_skycoord(ii,jj,WCS(hdr))
    print(point_coords)
    for axi in axs:
        axi.scatter(ii,jj,color='k',s=10)
    
    
    ax.text(0.0435,-300,r'$\ell$='+str(np.round(point_coords.l.deg,2))+r'$^{\circ}$'+
                       r' $b$='+str(np.round(point_coords.b.deg,2))+r'$^{\circ}$',fontsize=fs)
    ax.text(0.0435,-370,'GMIMS',fontsize=fs)
    ax.text(0.0435,-420,'CGPS',fontsize=fs)
    ax.text(0.0435,-470,'CGPS+GMIMS',fontsize=fs)
    
    ax.text(0.0441,-370,'RM='+str(np.round(RM_pt_G,1)),fontsize=fs)
    ax.text(0.0441,-420,'RM='+str(np.round(RM_pt_C,1)),fontsize=fs)
    ax.text(0.0441,-470,'RM='+str(np.round(RM_pt_CG,1)),fontsize=fs)
    
    ax.text(0.0446,-370,r'$\pm$ '+str(np.round(err_list[0][jj,ii],1))+r' rad m $^{-2}$',fontsize=fs)
    ax.text(0.0446,-420,r'$\pm$ '+str(np.round(err_list[1][jj,ii],1))+r' rad m $^{-2}$',fontsize=fs)
    ax.text(0.0446,-470,r'$\pm$ '+str(np.round(err_list[2][jj,ii],1))+r' rad m $^{-2}$',fontsize=fs)

    ax.text(0.0454,-370,r'$r=$'+str(np.round(r_list[0][jj,ii],1)),fontsize=fs)
    ax.text(0.0454,-420,r'$r=$'+str(np.round(r_list[1][jj,ii],1)),fontsize=fs)
    ax.text(0.0454,-470,r'$r=$'+str(np.round(r_list[2][jj,ii],1)),fontsize=fs)

    plt.show()

    return 
                                 
def make_the_plot(data_list,r_list,err_list,PAint_list,hdr,PA_G,PA_C,PA_CG,fs,
                  llim,blim,RMmax=100,PImax=0.5,*args,**kwargs):
  
    c = SkyCoord(llim, blim, frame=Galactic, unit="deg")
        
    fig1 = plt.figure(figsize=(20,10))
    
    ax0  = fig1.add_subplot(331, projection=WCS(hdr).celestial)
    im0  = ax0.imshow(data_list[0], origin='lower', vmin=-RMmax, vmax=RMmax,cmap='Spectral_r')
    ax0.set_xlim(WCS(hdr).world_to_pixel(c)[0])
    ax0.set_ylim(WCS(hdr).world_to_pixel(c)[1])
    #cbar1 = fig1.colorbar(im1,orientation='horizontal',fraction=0.06)
    
    ax1  = fig1.add_subplot(332, projection=WCS(hdr).celestial,sharex=ax0, sharey=ax0)
    im1  = ax1.imshow(data_list[1], origin='lower', vmin=-RMmax, vmax=RMmax,cmap='Spectral_r')
    
    ax2  = fig1.add_subplot(333, projection=WCS(hdr).celestial,sharex=ax0, sharey=ax0)
    im2  = ax2.imshow(data_list[2], origin='lower', vmin=0, vmax=PImax,cmap='cubehelix')
    
    ax3  = fig1.add_subplot(335, projection=WCS(hdr).celestial,sharex=ax0, sharey=ax0)
    im3  = ax3.imshow(data_list[3], origin='lower', vmin=-RMmax, vmax=RMmax,cmap='Spectral_r')
    
    ax4  = fig1.add_subplot(336, projection=WCS(hdr).celestial,sharex=ax0, sharey=ax0)
    im4  = ax4.imshow(data_list[4], origin='lower', vmin=0, vmax=PImax,cmap='cubehelix')
    
    ax5  = fig1.add_subplot(338, projection=WCS(hdr).celestial,sharex=ax0, sharey=ax0)
    im5  = ax5.imshow(data_list[5], origin='lower', vmin=-RMmax, vmax=RMmax,cmap='Spectral_r')
    #cbar5 = fig1.colorbar(im5,orientation='horizontal',fraction=0.08)
    
    ax6  = fig1.add_subplot(339, projection=WCS(hdr).celestial,sharex=ax0, sharey=ax0)
    im6  = ax6.imshow(data_list[6], origin='lower', vmin=0, vmax=PImax,cmap='cubehelix')
    #cbar6 = fig1.colorbar(im6,orientation='horizontal',fraction=0.08)
    
    for axs in [ax0,ax1,ax2,ax3,ax4,ax5,ax6]:
        axs.set_xlim(ax0.get_xlim())
        axs.set_ylim(ax0.get_ylim())
        axs.set_xlabel(' ')
        axs.set_ylabel(' ')
    
    #ax2  = fig1.add_subplot(223, projection=WCS(hdr).celestial, sharex=ax1, sharey=ax1)
    #im2  = ax2.imshow(data2, origin='lower', vmin=0, vmax=datamax2,cmap='cubehelix')
    #ax2.set_xlim(ax1.get_xlim())
    #ax2.set_ylim(ax1.get_ylim())
    #cbar2 = fig1.colorbar(im2,orientation='horizontal',fraction=0.06)
    #cbar2.set_label('PI (K)')
    
    ax  = fig1.add_subplot(334)
        
    cid = fig1.canvas.mpl_connect('button_press_event', 
                                  lambda event: mouse_event(event,ax,data_list,
                                                            r_list,err_list,PAint_list,hdr,
                                                            [ax0,ax1,ax2,ax3,ax4,ax5,ax6],
                                                            PA_G,PA_C,PA_CG,fs))
    klicker = clicker(ax5, ["s"], markers=[" "])
    ax5.get_legend().remove()
    
    plt.tight_layout()
    plt.show()
            
    return


if __name__ =='__main__': 
    do_RM_compare()