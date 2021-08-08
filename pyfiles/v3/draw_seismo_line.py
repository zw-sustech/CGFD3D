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



def draw_seismo_line(sac_fname_pre,sac_fname_end,recv_seq_no_start,recv_seq_no_end,recv_seq_no_incre,\
                     flag_show=1,figsize_line=[5,10],figdpi=150,flag_figsave=1,figpath='./fig',fignm_line='seismo_line.png'):

    # load data
    recv_indx = range(recv_seq_no_start,recv_seq_no_end+1,recv_seq_no_incre)
    nrecv = len(recv_indx)

    # receiver loop in line
    for irec in range(nrecv):

        sac_fname = sac_fname_pre + str(recv_indx[irec]) + sac_fname_end
    
        sac = read(sac_fname)
    
        sac_values = sac[0].data
        
        time_start = sac[0].stats['sac']['b']
        time_end   = sac[0].stats['sac']['e']
        time_npts  = sac[0].stats['npts']
        sac_time   = np.linspace(time_start,time_end,time_npts)
    
        if 'sac_line_values' not in locals():
            sac_line_values = np.empty((nrecv,time_npts))
        
        sac_line_values[irec,:] = sac_values

    
    # max amplitude for plotting scale
    scl_plt = np.max(np.abs(sac_line_values))
    
    # increment of show for Y-axis ticklabel
    ytick_incre = 9
    
    # create figure
    plt.figure(dpi = figdpi,figsize = (figsize_line[0],figsize_line[1]))
    
    # plot
    for irec in range(nrecv):
        plt.plot(sac_time,sac_line_values[irec,:]+irec*(2*scl_plt),'b',linewidth=1.0)
    
    plt.xlabel('Time (s)')
    plt.ylabel('Receiver ID')
    plt.title(sac_fname_pre[:-3])
    
    plt.xlim([np.min(sac_time),np.max(sac_time)])
    plt.ylim([-(1+0.5)*(2*scl_plt),(nrecv+0.5)*(2*scl_plt)])
    
    plt.yticks(np.arange(0,nrecv,ytick_incre)*(2*scl_plt),np.arange(1,nrecv+1,ytick_incre))
    
    # save to file
    if flag_figsave:
        subprocess.call('mkdir -p {}'.format(figpath),shell=True)
        figfullnm_line = figpath + '/' + fignm_line
        plt.savefig(figfullnm_line)
    
    # display on screen
    if flag_show:
        plt.show()


if __name__ == '__main__':

    # sac file name information
    sac_fname_pre = '../project/output/evt_test_single_force.line_x_1.no'
    sac_fname_end = '.Vz.sac'
    
    # receiver information to plot
    recv_seq_no_start = 0
    recv_seq_no_end   = 99
    recv_seq_no_incre = 2
    
    # figure control parameters
    # show figure on screen or not
    flag_show    = 1
    # figure size
    figsize_line = [5,10]
    # figure resolution
    figdpi       = 150
    # save figure or not
    flag_figsave = 1
    # figure path to save
    figpath      = './fig'
    # figure name to save
    fignm_line   = 'seismo_line.png'

    draw_seismo_line(sac_fname_pre,sac_fname_end,recv_seq_no_start,recv_seq_no_end,recv_seq_no_incre,\
                     flag_show,figsize_line,figdpi,flag_figsave,figpath,fignm_line)


