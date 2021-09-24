'''
Draw a cross-section of metric by using pcolor style
Author:      Yuanhang Huo
Email:       yhhuo@mail.ustc.edu.cn
Affiliation: University of Science and Technology of China
Date:        2021.06.24
'''

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys
sys.path.append(".")
from locate_metric import *
from gather_metric import *
from gather_coord import *

# ----------------------------- paremeters input -------------------------- #
# file and path name
#parfnm='./project/test.json'
#metric_dir='./project/output'
parfnm= '/home/zhangw/work/cgfd_arc/02src/test.json'
metric_dir='/home/zhangw/work/cgfd_arc/02src/output'

# which metric profile to plot
subs=[41,1,1]       # start from index '1'
subc=[1,-1,-1]      # '-1' to plot all points in this dimension
subt=[1,1,1]

# variable to plot
# 'jac', 'xi_x', 'xi_y', 'xi_z', 'eta_x', 'eta_y', 'eta_z',
# 'zeta_x', 'zeta_y', 'zeta_z'
varnm='jac'

# figure control parameters
# 1
flag_show    = 1
# 2
flag_figsave = 1
figpath      = './fig'
fignm        = 'metric.png'
figsize      = [4,4]
figdpi       = 300
# 3
flag_km     = 1
flag_title  = 1
clbtype     = 'jet'
clbrange    = [None,None]
# ------------------------------------------------------------------------- #

# load metric data
metricinfo=locate_metric(parfnm,'start',subs,'count',subc,'stride',subt,'metricdir',metric_dir)
# get coordinate data
[x,y,z]=gather_coord(metricinfo,'coorddir',metric_dir)
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

# load metric data
v=gather_metric(metricinfo,varnm,'metricdir',metric_dir)

# metric show
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
    plt.colorbar()

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
    plt.colorbar()

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
    plt.colorbar()


if flag_title:
    plt.title(varnm)

if flag_figsave:
    subprocess.call('mkdir -p {}'.format(figpath),shell=True)
    figfullnm=figpath + '/' + fignm
    plt.savefig(figfullnm)

if flag_show:
    plt.show()



