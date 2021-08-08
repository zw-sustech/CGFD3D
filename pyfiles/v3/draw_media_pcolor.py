'''
Draw a cross-section of media by using pcolor style
Author:      Yuanhang Huo
Email:       yhhuo@mail.ustc.edu.cn
Affiliation: University of Science and Technology of China
Date:        2021.08.08
'''

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys
sys.path.append(".")
from locate_media import *
from gather_media import *
from gather_coord import *

def draw_media_pcolor(parfnm,media_dir,subs,subc,subt,varnm,flag_show=1,flag_figsave=1,\
        figpath='./fig',fignm='media.png',figsize=[4,4],figdpi=150,flag_km=1,flag_title=1,\
        clbtype='jet',clbrange=[None,None]):

    # load media data
    mediainfo = locate_media(parfnm,'start',subs,'count',subc,'stride',subt,'mediadir',media_dir)
    # get coordinate data
    [x,y,z] = gather_coord(mediainfo,'coorddir',media_dir)
    nx = x.shape[0]
    ny = x.shape[1]
    nz = x.shape[2]
    # coordinate unit
    str_unit = 'm'
    if flag_km:
        x = x/1e3
        y = y/1e3
        z = z/1e3
        str_unit = 'km'
    
    # load media data
    if varnm == 'Vp':
        rho   = gather_media(mediainfo,'rho','mediadir',media_dir)
        mu    = gather_media(mediainfo,'mu','mediadir',media_dir)
        lambd = gather_media(mediainfo,'lambda','mediadir',media_dir)
        v = ((lambd+2*mu)/rho)**0.5
        v = v/1e3
    elif varnm == 'Vs':
        rho   = gather_media(mediainfo,'rho','mediadir',media_dir)
        mu    = gather_media(mediainfo,'mu','mediadir',media_dir)
        v = (mu/rho)**0.5
        v = v/1e3
    elif varnm == 'rho':
        v = gather_media(mediainfo,varnm,'mediadir',media_dir)
        v = v/1e3
    else:
        v = gather_media(mediainfo,varnm,'mediadir',media_dir)
    
    # media show
    plt.figure(dpi=figdpi,figsize=(figsize[0],figsize[1]))
    #plt.figure()
    
    if nx == 1:
    
        Y = np.squeeze(y).transpose(1,0)
        Y = np.row_stack((Y,Y[-1,:]))
        Y = np.column_stack((Y,Y[:,-1]+(Y[0,-1]-Y[0,-2])))
        Z = np.squeeze(z).transpose(1,0)
        Z = np.row_stack((Z,Z[-1,:]+(Z[-1,0]-Z[-2,0])))
        Z = np.column_stack((Z,Z[:,-1]))
        V = np.squeeze(v).transpose(1,0)
        V = np.row_stack((V,V[-1,:]))
        V = np.column_stack((V,V[:,-1]))
    
        plt.pcolor(Y,Z,V)
    
        plt.xlabel('Y ' + '(' + str_unit + ')')
        plt.ylabel('Z ' + '(' + str_unit + ')')
        plt.axis('image')
    
    elif ny == 1:
    
        X = np.squeeze(x).transpose(1,0)
        X = np.row_stack((X,X[-1,:]))
        X = np.column_stack((X,X[:,-1]+(X[0,-1]-X[0,-2])))
        Z = np.squeeze(z).transpose(1,0)
        Z = np.row_stack((Z,Z[-1,:]+(Z[-1,0]-Z[-2,0])))
        Z = np.column_stack((Z,Z[:,-1]))
        V = np.squeeze(v).transpose(1,0)
        V = np.row_stack((V,V[-1,:]))
        V = np.column_stack((V,V[:,-1]))
    
        plt.pcolor(X,Z,V)
    
        plt.xlabel('X ' + '(' + str_unit + ')')
        plt.ylabel('Z ' + '(' + str_unit + ')')
        plt.axis('image')
    
    else:
    
        X = np.squeeze(x).transpose(1,0)
        X = np.row_stack((X,X[-1,:]))
        X = np.column_stack((X,X[:,-1]+(X[0,-1]-X[0,-2])))
        Y = np.squeeze(y).transpose(1,0)
        Y = np.row_stack((Y,Y[-1,:]+(Y[-1,0]-Y[-2,0])))
        Y = np.column_stack((Y,Y[:,-1]))
        V = np.squeeze(v).transpose(1,0)
        V = np.row_stack((V,V[-1,:]))
        V = np.column_stack((V,V[:,-1]))
    
        plt.pcolor(X,Y,V)
    
        plt.xlabel('X ' + '(' + str_unit + ')')
        plt.ylabel('Y ' + '(' + str_unit + ')')
        plt.axis('image')
    
    
    clb = plt.colorbar()
    if varnm == 'Vp' or varnm == 'Vs':
        clb.ax.set_ylabel('$\mathrm{\mathsf{(km/s)}}$')
    if varnm == 'rho':
        clb.ax.set_ylabel('$\mathrm{\mathsf{(g/cm^3)}}$')
    
    if flag_title:
        plt.title(varnm)
    
    if flag_figsave:
        subprocess.call('mkdir -p {}'.format(figpath),shell=True)
        figfullnm=figpath + '/' + fignm
        plt.savefig(figfullnm)
    
    if flag_show:
        plt.show()


if __name__ == '__main__':
    
    # parameter json filename with path
    parfnm    = '../project/test.json'
    # media nc file path
    media_dir = '../project/output'
    
    # media starting index
    subs = [2,2,60]
    # media counting index
    subc = [-1,-1,1]
    # media stride index
    subt = [2,1,2]
    
    # variable name to plot
    varnm = 'Vp'
    
    # show figure or not
    flag_show    = 1
    # save figure or not
    flag_figsave = 1
    # figure path to save
    figpath = './fig'
    # figure name to save
    fignm   = 'media.png'
    # figure size to save
    figsize = [4,4]
    # figure resolution to save
    figdpi  = 150
    # axis unit km or m
    flag_km = 1
    # show figure title or not
    flag_title = 1
    # colorbar type
    clbtype    = 'jet'
    # colorbar range
    clbrange   = [None,None]

    draw_media_pcolor(parfnm,media_dir,subs,subc,subt,varnm,flag_show,flag_figsave,\
            figpath,fignm,figsize,figdpi,flag_km,flag_title,clbtype,clbrange)
    
