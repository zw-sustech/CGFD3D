'''
Draw a cross-section of media by using pcolor style
Author:      Yuanhang Huo
Email:       yhhuo@mail.ustc.edu.cn
Affiliation: University of Science and Technology of China
Date:        2021.06.23
'''

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(".")
from locate_media import *
from gather_media import *
from gather_coord import *

# ----------------------------- paremeters input -------------------------- #
# file and path name
parfnm='./project/test.json'
media_dir='./project/output'

# which media profile to plot
subs=[2,2,60]       # start from index '1'
subc=[-1,-1,1]      # '-1' to plot all points in this dimension
subt=[2,1,2]

# variable to plot
# 'Vp', 'Vs', 'rho', 'lambda', 'mu'
varnm='Vp'

# figure control parameters
# 1
flag_show    = 1
# 2
flag_figsave = 1
fignm        = 'media.png'
figsize      = [4,4]
figdpi       = 300
# 3
flag_km     = 1
flag_title  = 1
clbtype     = 'jet'
clbrange    = [None,None]
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
#plt.figure(dpi=figdpi,figsize=(figsize[0],figsize[1]))
plt.figure()

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
    plt.savefig(fignm)

if flag_show:
    plt.show()



