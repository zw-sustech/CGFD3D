'''
Plot the seismograms at a station
Author:      Yuanhang Huo
Email:       yhhuo@mail.ustc.edu.cn
Affiliation: University of Science and Technology of China
Date:        2021.06.25
'''

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import json
import sys
sys.path.append(".")
from obspy import read

# file and path name
sac_fname = '/home/zhangw/work/cgfd_arc/07/output/evt_by_par.pt_i50_j50_k60.Vx.sac'

# load data
sac = read(sac_fname)

sac_values = sac[0].data

time_start = sac[0].stats['sac']['b']
time_end   = sac[0].stats['sac']['e']
time_npts  = sac[0].stats['npts']
sac_time   = np.linspace(time_start, time_end, time_npts)

# create figure
figsize_rec  = [8,4]
figdpi       = 150
plt.figure(dpi=figdpi, figsize=(figsize_rec[0],figsize_rec[1]))

# plot
plt.plot(sac_time[:],sac_values[:],'b',linewidth=1.0)

plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title(sac_fname)

plt.xlim([np.min(sac_time[:]),np.max(sac_time[:])])

# save to file
flag_figsave = 1
figpath      = './fig'
fignm_rec    = 'seismo_rec.png'
if flag_figsave:
    subprocess.call('mkdir -p {}'.format(figpath),shell=True)
    figfullnm_rec=figpath + '/' + fignm_rec
    plt.savefig(figfullnm_rec)

# display on screen
flag_show = 1
if flag_show:
    plt.show()

