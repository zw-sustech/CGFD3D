'''
Draw a cross-section of media by using pcolor style
Author:      Yuanhang Huo
Email:       yhhuo@mail.ustc.edu.cn
Affiliation: University of Science and Technology of China
Date:        2021.06.23
'''

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import argparse
import sys
sys.path.append(".")
from locate_media import *
from gather_media import *
from gather_coord import *

# ---------------------------- parameters input --------------------------- #
parin=argparse.ArgumentParser(description='Introduction to the drawing script')
parin.add_argument('--parfnm',type=str,required=True,help='parameter json filename with path')
parin.add_argument('--media_dir',type=str,required=True,help='media nc file path')
parin.add_argument('--subs',type=str,required=True,help="starting index number (i,j,k) of media ('1' is the first index), e.g., [1,2,3]")
parin.add_argument('--subc',type=str,required=True,help='counting index number (i,j,k) of media, e.g., [10,1,30]')
parin.add_argument('--subt',type=str,required=True,help='stride number (i,j,k) of media, e.g., [2,2,2]')
parin.add_argument('--varnm',type=str,required=True,help="variable to plot, e.g., 'Vp','mu'")
parin.add_argument('--flag_show',type=int,default=1,help='show media or not, default=1')
parin.add_argument('--flag_figsave',type=int,default=1,help='save media figure or not, default=1')
parin.add_argument('--figpath',type=str,default='./fig',help='figure path to save, default=./fig')
parin.add_argument('--fignm',type=str,default='media.png',help='figure name to save, default=media.png')
parin.add_argument('--figsize',type=str,default='[4,4]',help='figure size to save, default=[4,4]')
parin.add_argument('--figdpi',type=int,default=300,help='figure resolution to save, default=300')
parin.add_argument('--flag_km',type=int,default=1,help='figure axis unit in km or not, default=1')
parin.add_argument('--flag_title',type=int,default=1,help='show figure title or not, default=1')
parin.add_argument('--clbtype',type=str,default='jet',help='colorbar type, default=jet')
parin.add_argument('--clbrange',type=str,default='[None,None]',help='colorbar range, default=[None,None]')

# get all input parameters
par=parin.parse_args()

# parameter json filename with path
parfnm=par.parfnm
# media nc file path
media_dir=par.media_dir

# media starting index
subsstr=par.subs.split(',')
subs=[int(subsstr[0][1:]),int(subsstr[1]),int(subsstr[2][:-1])]
# media counting index
subcstr=par.subc.split(',')
subc=[int(subcstr[0][1:]),int(subcstr[1]),int(subcstr[2][:-1])]
# media stride index
subtstr=par.subt.split(',')
subt=[int(subtstr[0][1:]),int(subtstr[1]),int(subtstr[2][:-1])]

# variable name to plot
varnm=par.varnm

# show figure or not
flag_show=par.flag_show
# save figure or not
flag_figsave=par.flag_figsave
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
# show figure title or not
flag_title=par.flag_title
# colorbar type
clbtype=par.clbtype
# colorbar range
clbrangestr=par.clbrange.split(',')
if clbrangestr[0][1:] == 'None' and clbrangestr[1][:-1] == 'None':
    clbrange=[None,None]
else:
    clbrange=[float(clbrangestr[0][1:]),float(clbrangestr[1][:-1])]

## print input parameters for QC
#print(parfnm,type(parfnm))
#print(media_dir,type(media_dir))
#print(subs,type(subs))
#print(subc,type(subc))
#print(subt,type(subt))
#print(varnm,type(varnm))
#print(flag_show,type(flag_show))
#print(flag_figsave,type(flag_figsave))
#print(figpath,type(figpath))
#print(fignm,type(fignm))
#print(figsize,type(figsize))
#print(figdpi,type(figdpi))
#print(flag_km,type(flag_km))
#print(flag_title,type(flag_title))
#print(clbtype,type(clbtype))
#print(clbrange,type(clbrange))
# ------------------------------------------------------------------------- #



# load media data
mediainfo=locate_media(parfnm,'start',subs,'count',subc,'stride',subt,'mediadir',media_dir)
# get coordinate data
[x,y,z]=gather_coord(mediainfo,'coorddir',media_dir)
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

# load media data
if varnm == 'Vp':
    rho=gather_media(mediainfo,'rho','mediadir',media_dir)
    mu=gather_media(mediainfo,'mu','mediadir',media_dir)
    lambd=gather_media(mediainfo,'lambda','mediadir',media_dir)
    v=((lambd+2*mu)/rho)**0.5
    v=v/1e3
elif varnm == 'Vs':
    rho=gather_media(mediainfo,'rho','mediadir',media_dir)
    mu=gather_media(mediainfo,'mu','mediadir',media_dir)
    v=(mu/rho)**0.5
    v=v/1e3
elif varnm == 'rho':
    v=gather_media(mediainfo,varnm,'mediadir',media_dir)
    v=v/1e3
else:
    v=gather_media(mediainfo,varnm,'mediadir',media_dir)

# media show
plt.figure(dpi=figdpi,figsize=(figsize[0],figsize[1]))
#plt.figure()

if nx == 1:

    Y=np.squeeze(y).transpose(1,0)
    Y=np.row_stack((Y,Y[-1,:]))
    Y=np.column_stack((Y,Y[:,-1]+(Y[0,-1]-Y[0,-2])))
    Z=np.squeeze(z).transpose(1,0)
    Z=np.row_stack((Z,Z[-1,:]+(Z[-1,0]-Z[-2,0])))
    Z=np.column_stack((Z,Z[:,-1]))
    V=np.squeeze(v).transpose(1,0)
    V=np.row_stack((V,V[-1,:]))
    V=np.column_stack((V,V[:,-1]))

    plt.pcolor(Y,Z,V)

    plt.xlabel('Y ' + '(' + str_unit + ')')
    plt.ylabel('Z ' + '(' + str_unit + ')')
    plt.axis('image')

elif ny == 1:

    X=np.squeeze(x).transpose(1,0)
    X=np.row_stack((X,X[-1,:]))
    X=np.column_stack((X,X[:,-1]+(X[0,-1]-X[0,-2])))
    Z=np.squeeze(z).transpose(1,0)
    Z=np.row_stack((Z,Z[-1,:]+(Z[-1,0]-Z[-2,0])))
    Z=np.column_stack((Z,Z[:,-1]))
    V=np.squeeze(v).transpose(1,0)
    V=np.row_stack((V,V[-1,:]))
    V=np.column_stack((V,V[:,-1]))

    plt.pcolor(X,Z,V)

    plt.xlabel('X ' + '(' + str_unit + ')')
    plt.ylabel('Z ' + '(' + str_unit + ')')
    plt.axis('image')

else:

    X=np.squeeze(x).transpose(1,0)
    X=np.row_stack((X,X[-1,:]))
    X=np.column_stack((X,X[:,-1]+(X[0,-1]-X[0,-2])))
    Y=np.squeeze(y).transpose(1,0)
    Y=np.row_stack((Y,Y[-1,:]+(Y[-1,0]-Y[-2,0])))
    Y=np.column_stack((Y,Y[:,-1]))
    V=np.squeeze(v).transpose(1,0)
    V=np.row_stack((V,V[-1,:]))
    V=np.column_stack((V,V[:,-1]))

    plt.pcolor(X,Y,V)

    plt.xlabel('X ' + '(' + str_unit + ')')
    plt.ylabel('Y ' + '(' + str_unit + ')')
    plt.axis('image')


clb=plt.colorbar()
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



