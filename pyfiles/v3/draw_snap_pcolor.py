'''
Draw seismic wavefield snapshot by using pcolor style
Author:      Yuanhang Huo
Email:       yhhuo@mail.ustc.edu.cn
Affiliation: University of Science and Technology of China
Date:        2021.08.08
'''

import numpy as np
import matplotlib.pyplot as plt
import imageio
import subprocess
import sys
sys.path.append(".")
from locate_snap import *
from gather_snap import *
from gather_coord import *

def draw_snap_pcolor(parfnm,coord_dir,snap_dir,id,subs,subc,subt,varnm,ns,ne,nt,\
        flag_show=1,taut=0.5,flag_imgsave=1,flag_gifsave=1,figpath='./fig',\
        fignm='fd3dsnap.png',figsize=[4,4],figdpi=150,flag_km=1,\
        clbtype='seismic',clbrange=[None,None]):

    # load snapshot data
    snapinfo = locate_snap(parfnm,id,'start',subs,'count',subc,'stride',subt,'snapdir',snap_dir)
    
    # get coordinate data
    [x,y,z] = gather_coord(snapinfo,'coorddir',coord_dir)
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
    
    # snapshot show
    plt.figure(dpi=figdpi,figsize=(figsize[0],figsize[1]))
    #plt.figure()
    
    for nlayer in range(ns,ne+nt,nt):
    
        v,t = gather_snap(snapinfo,nlayer,varnm,'snapdir',snap_dir)
        print("Drawing " + str(nlayer) + "th time step (t = " + "{:.4f}".format(t) + " s)")
    
        figtitle = "Snapshot of " + varnm + " at " + "{:9.4f}".format(t) + " s"
    
        if nx == 1:
    
            #Y = np.squeeze(y).transpose(1,0)
            #Y = np.row_stack((Y,Y[-1,:]))
            #Y = np.column_stack((Y,Y[:,-1]+(Y[0,-1]-Y[0,-2])))
            #Z = np.squeeze(z).transpose(1,0)
            #Z = np.row_stack((Z,Z[-1,:]+(Z[-1,0]-Z[-2,0])))
            #Z = np.column_stack((Z,Z[:,-1]))
            #V = np.squeeze(v).transpose(1,0)
            #V = np.row_stack((V,V[-1,:]))
            #V = np.column_stack((V,V[:,-1]))
            
            Y = np.squeeze(y).transpose(1,0)
            Z = np.squeeze(z).transpose(1,0)
            V = np.squeeze(v).transpose(1,0)
    
            plt.clf()
            plt.cla()
            plt.pcolor(Y,Z,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
    
            plt.xlabel('Y ' + '(' + str_unit + ')')
            plt.ylabel('Z ' + '(' + str_unit + ')')
            plt.title(figtitle)
            plt.colorbar()
            plt.axis('image')
            if flag_show:
                plt.pause(taut)
    
        elif ny == 1:
            
            #X = np.squeeze(x).transpose(1,0)
            #X = np.row_stack((X,X[-1,:]))
            #X = np.column_stack((X,X[:,-1]+(X[0,-1]-X[0,-2])))
            #Z = np.squeeze(z).transpose(1,0)
            #Z = np.row_stack((Z,Z[-1,:]+(Z[-1,0]-Z[-2,0])))
            #Z = np.column_stack((Z,Z[:,-1]))
            #V = np.squeeze(v).transpose(1,0)
            #V = np.row_stack((V,V[-1,:]))
            #V = np.column_stack((V,V[:,-1]))
    
            X = np.squeeze(x).transpose(1,0)
            Z = np.squeeze(z).transpose(1,0)
            V = np.squeeze(v).transpose(1,0)
    
            plt.clf()
            plt.cla()
            plt.pcolor(X,Z,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
    
            plt.xlabel('X ' + '(' + str_unit + ')')
            plt.ylabel('Z ' + '(' + str_unit + ')')
            plt.title(figtitle)
            plt.colorbar()
            plt.axis('image')
            if flag_show:
                plt.pause(taut)
    
        else:
    
            #X = np.squeeze(x).transpose(1,0)
            #X = np.row_stack((X,X[-1,:]))
            #X = np.column_stack((X,X[:,-1]+(X[0,-1]-X[0,-2])))
            #Y = np.squeeze(y).transpose(1,0)
            #Y = np.row_stack((Y,Y[-1,:]+(Y[-1,0]-Y[-2,0])))
            #Y = np.column_stack((Y,Y[:,-1]))
            #V = np.squeeze(v).transpose(1,0)
            #V = np.row_stack((V,V[-1,:]))
            #V = np.column_stack((V,V[:,-1]))
    
            X = np.squeeze(x).transpose(1,0)
            Y = np.squeeze(y).transpose(1,0)
            V = np.squeeze(v).transpose(1,0)
    
            plt.clf()
            plt.cla()
            plt.pcolor(X,Y,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])
    
            plt.xlabel('X ' + '(' + str_unit + ')')
            plt.ylabel('Y ' + '(' + str_unit + ')')
            plt.title(figtitle)
            plt.colorbar()
            plt.axis('image')
            if flag_show:
                plt.pause(taut)
    
        if flag_imgsave or flag_gifsave:
            subprocess.call('mkdir -p {}'.format(figpath),shell=True)
            imgnm  = fignm[:-(len(fignm.split('.')[-1])+1)]
            imgfmt = fignm.split('.')[-1]
            imgfullnm = '{}/{}_timestep_{}.{}'.format(figpath,imgnm,nlayer,imgfmt)
            plt.savefig(imgfullnm)
    
    if flag_gifsave:
        frames = []
        for nlayer in range(ns,ne+nt,nt):
            imgfullnm = '{}/{}_timestep_{}.{}'.format(figpath,imgnm,nlayer,imgfmt)
            frames.append(imageio.imread(imgfullnm))
        imageio.mimsave('{}/{}.gif'.format(figpath,imgnm),frames,'GIF',duration=taut)
    
    if flag_imgsave == 0 and flag_gifsave == 1:
        subprocess.call('rm {}/{}_timestep_*.{}'.format(figpath,imgnm,imgfmt),shell=True)
    
    if flag_show:
        plt.show()


if __name__ == '__main__':
    
    # parameter json filename with path
    parfnm = '/home/zhangw/work/wpsfd_dt/02/test.json'
    # snapshot nc file path
    snap_dir = '/home/zhangw/work/wpsfd_dt/02/output'
    # coord nc file path
    coord_dir = '/home/zhangw/work/wpsfd_dt/02/output'
    
    # snapshot id
    id = 1
    ## snapshot starting index
    #subs = [1,1,50]
    ## snapshot counting index
    #subc = [-1,-1,1]
    ## snapshot stride index
    #subt = [1,1,1]

    #-- ix slice
    #subs = [41,1,1]
    #subc = [1,-1,-1]
    #subt = [1,1,1]
    #figpath  = './fig.ix'

    #-- iy slice
    #subs = [1,41,1]
    #subc = [-1,1,-1]
    #subt = [1,1,1]
    #figpath  = './fig.iy5k'

    #-- iz slice
    subs = [1,1,51]
    subc = [-1,-1,1]
    subt = [1,1,1]
    figpath  = './fig.iz51'
    
    # variable name to plot
    varnm = 'Vz'
    # starting time step to plot
    ns = 1
    # ending time step to plot
    ne = 500
    # time stride to plot
    nt = 10
    
    # show figure or not
    flag_show = 1
    # plotting pause time in second
    taut = 0.5
    # save figure or not
    flag_imgsave = 1
    # save gif or not
    flag_gifsave = 0
    # figure name to save
    fignm    = 'fd3dsnap.png'
    # figure size to save
    figsize  = [4,4]
    # figure resolution to save
    figdpi   = 150
    # axis unit km or m
    flag_km  = 1
    # colorbar type
    clbtype  = 'seismic'
    # colorbar range
    clbrange = [None,None]

    draw_snap_pcolor(parfnm,coord_dir,snap_dir,id,subs,subc,subt,varnm,ns,ne,nt,\
            flag_show,taut,flag_imgsave,flag_gifsave,figpath,fignm,figsize,\
            figdpi,flag_km,clbtype,clbrange)
