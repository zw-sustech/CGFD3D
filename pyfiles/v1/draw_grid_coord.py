'''
Plot the grid
Author:      Yuanhang Huo
Email:       yhhuo@mail.ustc.edu.cn
Affiliation: University of Science and Technology of China
Date:        2021.06.23
'''

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys
sys.path.append(".")
from locate_coord import *
from gather_coord import *

# ---------------------------- parameters input --------------------------- #
# file and path name
#parfnm='./project/test.json'
#coord_dir='./project/output'
parfnm= '/home/zhangw/work/cgfd_arc/02src/test.json'
coord_dir='/home/zhangw/work/cgfd_arc/02src/output'

# which grid profile to plot
subs=[1,1,50]        # start from index '1'
subc=[-1,-1,1]      # '-1' to plot all points in this dimension
subt=[1,1,1]

# x increment to plot
pltincre_x = 2
# y increment to plot
pltincre_y = 2

# figure control parameters
# 1
flag_show    = 1
# 2
flag_figsave = 1
figpath      = './fig'
fignm        = 'grid.png'
figsize      = [4,4]
figdpi       = 300
# 3
flag_km     = 1
flag_title  = 1
# ------------------------------------------------------------------------- #



# load grid coordinate
coordinfo=locate_coord(parfnm,'start',subs,'count',subc,'stride',subt,'coorddir',coord_dir)
[x,y,z]=gather_coord(coordinfo,'coorddir',coord_dir)
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

# plot
x=np.squeeze(x)
y=np.squeeze(y)
z=np.squeeze(z)

# grid show
plt.figure(dpi=figdpi,figsize=(figsize[0],figsize[1]))
#plt.figure()

if nx == 1:

    plt.plot(y[::pltincre_x,::pltincre_y].transpose(1,0),\
             z[::pltincre_x,::pltincre_y].transpose(1,0),\
             'k-')
    plt.plot(y[::pltincre_x,::pltincre_y],\
             z[::pltincre_x,::pltincre_y],\
             'k-')

    plt.xlabel('Y axis (' + str_unit + ')')
    plt.xlabel('Z axis (' + str_unit + ')')
    plt.axis('image')

elif ny == 1:

    plt.plot(x[::pltincre_x,::pltincre_y].transpose(1,0),\
             z[::pltincre_x,::pltincre_y].transpose(1,0),\
             'k-')
    plt.plot(x[::pltincre_x,::pltincre_y],\
             z[::pltincre_x,::pltincre_y],\
             'k-')

    plt.xlabel('X axis (' + str_unit + ')')
    plt.xlabel('Z axis (' + str_unit + ')')
    plt.axis('image')

else:

    plt.plot(x[::pltincre_x,::pltincre_y].transpose(1,0),\
             y[::pltincre_x,::pltincre_y].transpose(1,0),\
             'k-')
    plt.plot(x[::pltincre_x,::pltincre_y],\
             y[::pltincre_x,::pltincre_y],\
             'k-')

    plt.xlabel('X axis (' + str_unit + ')')
    plt.xlabel('Y axis (' + str_unit + ')')
    plt.axis('image')

# title
if flag_title:
    if nx == 1:
        gridtitle='YOZ-Grid'
    elif ny == 1:
        gridtitle='XOZ-Grid'
    else:
        gridtitle='XOY-Grid'
    plt.title(gridtitle)

if flag_figsave:
    subprocess.call('mkdir -p {}'.format(figpath),shell=True)
    figfullnm=figpath + '/' + fignm
    plt.savefig(figfullnm)

if flag_show:
    plt.show()



