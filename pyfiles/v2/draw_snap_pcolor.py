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
import argparse
import sys
sys.path.append(".")
from locate_snap import *
from gather_snap import *
from gather_coord import *

# ---------------------------- parameters input --------------------------- #
parin=argparse.ArgumentParser(description='Introduction to the drawing script')
parin.add_argument('--parfnm',type=str,required=True,help='parameter json filename with path')
parin.add_argument('--snap_dir',type=str,required=True,help='snapshot nc file path')
parin.add_argument('--id',type=int,required=True,help='snapshot id')
parin.add_argument('--subs',type=str,required=True,help="starting index number (i,j,k) of snapshot ('1' is the first index), e.g., [1,2,3]")
parin.add_argument('--subc',type=str,required=True,help='counting index number (i,j,k) of snapshot, e.g., [10,1,30]')
parin.add_argument('--subt',type=str,required=True,help='stride number (i,j,k) of snapshot, e.g., [2,2,2]')
parin.add_argument('--varnm',type=str,required=True,help="variable to plot, e.g., 'Vx','Tzz'")
parin.add_argument('--ns',type=int,required=True,help='starting time step to plot')
parin.add_argument('--ne',type=int,required=True,help='ending time step to plot')
parin.add_argument('--nt',type=int,required=True,help='time stride to plot')
parin.add_argument('--flag_show',type=int,default=1,help='show snapshot or not, default=1')
parin.add_argument('--taut',type=float,default=0.5,help='plotting pause time (in second) between two time steps, default=0.5')
parin.add_argument('--flag_imgsave',type=int,default=1,help='save snapshot figure or not, default=1')
parin.add_argument('--flag_gifsave',type=int,default=1,help='save snapshot gif or not, default=1')
parin.add_argument('--figpath',type=str,default='./fig',help='figure path to save, default=./fig')
parin.add_argument('--fignm',type=str,default='fd3dsnap.png',help='figure name to save, default=fd3dsnap.png')
parin.add_argument('--figsize',type=str,default='[4,4]',help='figure size to save, default=[4,4]')
parin.add_argument('--figdpi',type=int,default=300,help='figure resolution to save, default=300')
parin.add_argument('--flag_km',type=int,default=1,help='figure axis unit in km or not, default=1')
parin.add_argument('--clbtype',type=str,default='seismic',help='colorbar type, default=seismic')
parin.add_argument('--clbrange',type=str,default='[None,None]',help='colorbar range, default=[None,None]')

# get all input parameters
par=parin.parse_args()

# parameter json filename with path
parfnm=par.parfnm
# snapshot nc file path
snap_dir=par.snap_dir

# snapshot id
id=par.id
# snapshot starting index
subsstr=par.subs.split(',')
subs=[int(subsstr[0][1:]),int(subsstr[1]),int(subsstr[2][:-1])]
# snapshot counting index
subcstr=par.subc.split(',')
subc=[int(subcstr[0][1:]),int(subcstr[1]),int(subcstr[2][:-1])]
# snapshot stride index
subtstr=par.subt.split(',')
subt=[int(subtstr[0][1:]),int(subtstr[1]),int(subtstr[2][:-1])]

# variable name to plot
varnm=par.varnm
# starting time step to plot
ns=par.ns
# ending time step to plot
ne=par.ne
# time stride to plot
nt=par.nt

# show figure or not
flag_show=par.flag_show
# plotting pause time in second
taut=par.taut
# save figure or not
flag_imgsave=par.flag_imgsave
# save gif or not
flag_gifsave=par.flag_gifsave
# figure path to save
figpath=par.figpath
# figure name to save
fignm=par.fignm
# figure size to save
figsizestr=par.figsize.split(',')
figsize=[int(figsizestr[0][1:]),int(figsizestr[1][:-1])]
# figure resolution to save
figdpi=par.figdpi
# axis unit km or m
flag_km=par.flag_km
# colorbar type
clbtype=par.clbtype
# colorbar range
clbrangestr=par.clbrange.split(',')
if clbrangestr[0][1:] == 'None' and clbrangestr[1][:-1] == 'None':
    clbrange=[None,None]
else:
    clbrange=[float(clbrangestr[0][1:]),float(clbrangestr[1][:-1])]

# print input parameters for QC
#print(parfnm,type(parfnm))
#print(snap_dir,type(snap_dir))
#print(id,type(id))
#print(subs,type(subs))
#print(subc,type(subc))
#print(subt,type(subt))
#print(varnm,type(varnm))
#print(ns,type(ns))
#print(ne,type(ne))
#print(nt,type(nt))
#print(flag_show,type(flag_show))
#print(taut,type(taut))
#print(flag_imgsave,type(flag_imgsave))
#print(flag_gifsave,type(flag_gifsave))
#print(figpath,type(figpath))
#print(fignm,type(fignm))
#print(figsize,type(figsize))
#print(figdpi,type(figdpi))
#print(flag_km,type(flag_km))
#print(clbtype,type(clbtype))
#print(clbrange,type(clbrange))
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
plt.figure(dpi=figdpi,figsize=(figsize[0],figsize[1]))
#plt.figure()

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

