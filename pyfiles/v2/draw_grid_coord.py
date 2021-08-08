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
import argparse
import sys
sys.path.append(".")
from locate_coord import *
from gather_coord import *

# ---------------------------- parameters input --------------------------- #
parin=argparse.ArgumentParser(description='Introduction to the drawing script')
parin.add_argument('--parfnm',type=str,required=True,help='parameter json filename with path')
parin.add_argument('--coord_dir',type=str,required=True,help='coordinate nc file path')
parin.add_argument('--subs',type=str,required=True,help="starting index number (i,j,k) of coordinate ('1' is the first index), e.g., [1,2,3]")
parin.add_argument('--subc',type=str,required=True,help='counting index number (i,j,k) of coordinate, e.g., [10,1,30]')
parin.add_argument('--subt',type=str,required=True,help='stride number (i,j,k) of coordinate, e.g., [2,2,2]')
parin.add_argument('--pltincre_x',type=int,default=2,help='x increment to plot, default=2')
parin.add_argument('--pltincre_y',type=int,default=2,help='y increment to plot, default=2')
parin.add_argument('--flag_show',type=int,default=1,help='show figure or not, default=1')
parin.add_argument('--flag_figsave',type=int,default=1,help='save grid figure or not, default=1')
parin.add_argument('--figpath',type=str,default='./fig',help='figure path to save, default=./fig')
parin.add_argument('--fignm',type=str,default='grid.png',help='figure name to save, default=grid.png')
parin.add_argument('--figsize',type=str,default='[4,4]',help='figure size to save, default=[4,4]')
parin.add_argument('--figdpi',type=int,default=300,help='figure resolution to save, default=300')
parin.add_argument('--flag_km',type=int,default=1,help='figure axis unit in km or not, default=1')
parin.add_argument('--flag_title',type=int,default=1,help='show figure title or not, default=1')

# get all input parameters
par=parin.parse_args()

# parameter json filename with path
parfnm=par.parfnm
# coordinate nc file path
coord_dir=par.coord_dir

# coordinate starting index
subsstr=par.subs.split(',')
subs=[int(subsstr[0][1:]),int(subsstr[1]),int(subsstr[2][:-1])]
# coordinate counting index
subcstr=par.subc.split(',')
subc=[int(subcstr[0][1:]),int(subcstr[1]),int(subcstr[2][:-1])]
# coordinate stride index
subtstr=par.subt.split(',')
subt=[int(subtstr[0][1:]),int(subtstr[1]),int(subtstr[2][:-1])]

# x increment to plot
pltincre_x=par.pltincre_x
# y increment to plot
pltincre_y=par.pltincre_y

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
# show title or not
flag_title=par.flag_title

## print input parameters for QC
#print(parfnm,type(parfnm))
#print(coord_dir,type(coord_dir))
#print(subs,type(subs))
#print(subc,type(subc))
#print(subt,type(subt))
#print(pltincre_x,type(pltincre_x))
#print(pltincre_y,type(pltincre_y))
#print(flag_show,type(flag_show))
#print(flag_figsave,type(flag_figsave))
#print(figpath,type(figpath))
#print(fignm,type(fignm))
#print(figsize,type(figsize))
#print(figdpi,type(figdpi))
#print(flag_km,type(flag_km))
#print(flag_title,type(flag_title))
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



