'''
Draw seismic wavefield snapshot by using pcolor style
Author:      Yuanhang Huo
Email:       yhhuo@mail.ustc.edu.cn
Affiliation: University of Science and Technology of China
Date:        2021.06.23
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

# ---------------------------- parameters input --------------------------- #
# file and path name
parfnm= '/home/zhangw/work/cgfd_arc/09/test.json'
snap_dir='/home/zhangw/work/cgfd_arc/09/output'

# which snapshot to plot
id=1

#- z plane
#subs=[1,1,50]       # start from index '1'
#subc=[-1,-1,1]      # '-1' to plot all points in this dimension
#subt=[1,1,1]

#- x plane
subs=[41,1,1]       # start from index '1'
subc=[1,-1,-1]      # '-1' to plot all points in this dimension
subt=[1,1,1]

#- y plane
#subs=[1,41,1]       # start from index '1'
#subc=[-1,1,-1]      # '-1' to plot all points in this dimension
#subt=[1,1,1]

# variable and time to plot
varnm='Vz'
ns=1
ne=500
nt=10

# figure control parameters
# 1
flag_show    = 1
taut         = 0.5
# 2
flag_imgsave = 1
flag_gifsave = 1
figpath      = './fig'
fignm        = 'fd3dsnap.png'
figsize      = [4,4]
figdpi       = 150
# 3
flag_km      = 1
clbtype      = 'seismic'
clbrange     = [None,None]
# ------------------------------------------------------------------------- #



# load snapshot data
snapinfo=locate_snap(parfnm,id,'start',subs,'count',subc,'stride',subt,'snapdir',snap_dir)

# get coordinate data
[x,y,z]=gather_coord(snapinfo,'coorddir',snap_dir)
nx=x.shape[0]
ny=x.shape[1]
nz=x.shape[2]
# coordinate unit
str_unit='m'
if flag_km:
    x=x/1e3
    y=y/1e3
    z=z/1e3
    str_unit='km'

# snapshot show
#plt.figure(dpi=figdpi,figsize=(figsize[0],figsize[1]))
plt.figure()

for nlayer in range(ns,ne+nt,nt):

    v,t=gather_snap(snapinfo,nlayer,varnm,'snapdir',snap_dir)
    print("Drawing " + str(nlayer) + "th time step (t = " + "{:.4f}".format(t) + " s)")

    figtitle="Snapshot of " + varnm + " at " + "{:9.4f}".format(t) + " s"

    if nx == 1:

        #Y=np.squeeze(y).transpose(1,0)
        #Y=np.row_stack((Y,Y[-1,:]))
        #Y=np.column_stack((Y,Y[:,-1]+(Y[0,-1]-Y[0,-2])))
        #Z=np.squeeze(z).transpose(1,0)
        #Z=np.row_stack((Z,Z[-1,:]+(Z[-1,0]-Z[-2,0])))
        #Z=np.column_stack((Z,Z[:,-1]))
        #V=np.squeeze(v).transpose(1,0)
        #V=np.row_stack((V,V[-1,:]))
        #V=np.column_stack((V,V[:,-1]))
        
        Y=np.squeeze(y).transpose(1,0)
        Z=np.squeeze(z).transpose(1,0)
        V=np.squeeze(v).transpose(1,0)

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
        
        #X=np.squeeze(x).transpose(1,0)
        #X=np.row_stack((X,X[-1,:]))
        #X=np.column_stack((X,X[:,-1]+(X[0,-1]-X[0,-2])))
        #Z=np.squeeze(z).transpose(1,0)
        #Z=np.row_stack((Z,Z[-1,:]+(Z[-1,0]-Z[-2,0])))
        #Z=np.column_stack((Z,Z[:,-1]))
        #V=np.squeeze(v).transpose(1,0)
        #V=np.row_stack((V,V[-1,:]))
        #V=np.column_stack((V,V[:,-1]))

        X=np.squeeze(x).transpose(1,0)
        Z=np.squeeze(z).transpose(1,0)
        V=np.squeeze(v).transpose(1,0)

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

        #X=np.squeeze(x).transpose(1,0)
        #X=np.row_stack((X,X[-1,:]))
        #X=np.column_stack((X,X[:,-1]+(X[0,-1]-X[0,-2])))
        #Y=np.squeeze(y).transpose(1,0)
        #Y=np.row_stack((Y,Y[-1,:]+(Y[-1,0]-Y[-2,0])))
        #Y=np.column_stack((Y,Y[:,-1]))
        #V=np.squeeze(v).transpose(1,0)
        #V=np.row_stack((V,V[-1,:]))
        #V=np.column_stack((V,V[:,-1]))

        X=np.squeeze(x).transpose(1,0)
        Y=np.squeeze(y).transpose(1,0)
        V=np.squeeze(v).transpose(1,0)

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
        imgnm=fignm[:-(len(fignm.split('.')[-1])+1)]
        imgfmt=fignm.split('.')[-1]
        imgfullnm='{}/{}_timestep_{}.{}'.format(figpath,imgnm,nlayer,imgfmt)
        plt.savefig(imgfullnm)

if flag_gifsave:
    frames=[]
    for nlayer in range(ns,ne+nt,nt):
        imgfullnm='{}/{}_timestep_{}.{}'.format(figpath,imgnm,nlayer,imgfmt)
        frames.append(imageio.imread(imgfullnm))
    imageio.mimsave('{}/{}.gif'.format(figpath,imgnm),frames,'GIF',duration=taut)

if flag_imgsave == 0 and flag_gifsave == 1:
    subprocess.call('rm {}/{}_timestep_*.{}'.format(figpath,imgnm,imgfmt),shell=True)
        

if flag_show:
    plt.show()


