'''
Plot the seismograms at a station
Author:      Yuanhang Huo
Email:       yhhuo@mail.ustc.edu.cn
Affiliation: University of Science and Technology of China
Date:        2021.08.07
'''

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import json
import sys
sys.path.append(".")
from obspy import read


def draw_seismo_station(sac_fname,flag_show=1,figsize_rec=[8,4],figdpi=150,flag_figsave=1,figpath='./fig',fignm_rec='seismo_rec.png'):

    # load data
    sac = read(sac_fname)
    
    sac_values = sac[0].data
    
    time_start = sac[0].stats['sac']['b']
    time_end   = sac[0].stats['sac']['e']
    time_npts  = sac[0].stats['npts']
    sac_time   = np.linspace(time_start, time_end, time_npts)
    
    # create figure
    plt.figure(dpi=figdpi, figsize=(figsize_rec[0],figsize_rec[1]))
    
    # plot
    plt.plot(sac_time[:],sac_values[:],'b',linewidth=1.0)
    
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title(sac_fname)
    
    plt.xlim([np.min(sac_time[:]),np.max(sac_time[:])])
    
    # save to file
    if flag_figsave:
        subprocess.call('mkdir -p {}'.format(figpath),shell=True)
        figfullnm_rec = figpath + '/' + fignm_rec
        plt.savefig(figfullnm_rec)
    
    # display on screen
    if flag_show:
        plt.show()


if __name__ == '__main__':

    # file and path name
    sac_fname = '../project/output/evt_test_single_force.line_x_1.no68.Vz.sac'
    
    # figure control parameters
    # show figure on screen or not
    flag_show = 1
    # figure size
    figsize_rec  = [8,4]
    # figure resolution
    figdpi       = 150
    # save figure or not
    flag_figsave = 1
    # figure path to save
    figpath      = './fig'
    # figure name to save
    fignm_rec    = 'seismo_rec.png'

    draw_seismo_station(sac_fname,flag_show,figsize_rec,figdpi,flag_figsave,figpath,fignm_rec)



