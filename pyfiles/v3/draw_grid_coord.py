'''
Plot the grid
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
from locate_coord import *
from gather_coord import *

def draw_grid_coord(parfnm,coord_dir,subs,subc,subt,pltincre_x=2,pltincre_y=2,\
        flag_show=1,flag_figsave=1,figpath='./fig',fignm='grid.png',figsize=[4,4],\
        figdpi=150,flag_km=1,flag_title=1):

    # load grid coordinate
    coordinfo = locate_coord(parfnm,'start',subs,'count',subc,'stride',subt,'coorddir',coord_dir)
    [x,y,z] = gather_coord(coordinfo,'coorddir',coord_dir)
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
    
    # plot
    x = np.squeeze(x)
    y = np.squeeze(y)
    z = np.squeeze(z)
    
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
            gridtitle = 'YOZ-Grid'
        elif ny == 1:
            gridtitle = 'XOZ-Grid'
        else:
            gridtitle = 'XOY-Grid'
        plt.title(gridtitle)
    
    if flag_figsave:
        subprocess.call('mkdir -p {}'.format(figpath),shell=True)
        figfullnm=figpath + '/' + fignm
        plt.savefig(figfullnm)
    
    if flag_show:
        plt.show()


if __name__ == '__main__':
    
    # parameter json filename with path
    parfnm    = '../project/test.json'
    # coordinate nc file path
    coord_dir = '../project/output'
    
    # coordinate starting index
    subs = [5,5,1]
    # coordinate counting index
    subc = [-1,-1,1]
    # coordinate stride index
    subt = [1,1,1]

    # x increment to plot
    pltincre_x = 2
    # y increment to plot
    pltincre_y = 2
    
    # show figure or not
    flag_show    = 1
    # save figure or not
    flag_figsave = 1
    # figure path to save
    figpath    = './fig'
    # figure name to save
    fignm      = 'grid.png'
    # figure size to save
    figsize    = [4,4]
    # figure resolution to save
    figdpi     = 150
    # axis unit km or m
    flag_km    = 1
    # show title or not
    flag_title = 1

    draw_grid_coord(parfnm,coord_dir,subs,subc,subt,pltincre_x,pltincre_y,\
            flag_show,flag_figsave,figpath,fignm,figsize,figdpi,flag_km,flag_title)

