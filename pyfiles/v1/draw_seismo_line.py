'''
Plot the seismograms at a line
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

# --------------------------- parameters input ---------------------------- #
# file and path name
parfnm='./project/test.json'
seismo_dir='./project/output'

# which line to plot (start from index '1')
lineid=1

# which receiver of the line to plot (start from index '0')
recid=0

# which variable to plot
varnm='Vz'

# figure control paramters
# 1
flag_show    = 1
# 2
flag_figsave = 1
figpath      = './fig'
fignm_rec    = 'seismo_rec.png'
figsize_rec  = [8,4]
fignm_line   = 'seismo_line.png'
figsize_line = [5,10]
figdpi       = 150
# ------------------------------------------------------------------------- #


# raed parameter file
par=json.load(open(parfnm,'r'))

# judge source type
if 'single_force' in par['source_input'].keys():
    lineprefix=par['source_input']['single_force']['name']
elif 'single_moment' in par['source_input'].keys():
    lineprefix=par['source_input']['single_moment']['name']
else:
    lineprefix='user_input_src'

# line name and receiver number
linenm=par['receiver_line'][lineid-1]['name']
nrec=par['receiver_line'][lineid-1]['grid_index_count']

# load data
for irec in range(nrec):
    sacnm=seismo_dir + '/' + lineprefix + '.' + linenm + '.' + 'no' + str(irec) + '.' + varnm + '.sac'
    sacfile=read(sacnm)
    if 'seismodata' not in locals():
        seismodata=np.empty((nrec,sacfile[0].stats['npts']))
    if 'seismot' not in locals():
        seismot=np.empty((nrec,sacfile[0].stats['npts']))
    seismodata[irec,:]=sacfile[0].data
    seismot[irec,:]=np.linspace(sacfile[0].stats['sac']['b'],sacfile[0].stats['sac']['e'],\
                                sacfile[0].stats['npts'])

# plot single receiver
#plt.figure()
plt.figure(dpi=figdpi,figsize=(figsize_rec[0],figsize_rec[1]))
plt.plot(seismot[recid,:],seismodata[recid,:],'b',linewidth=1.0)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title(varnm + ' at No.' + str(recid+1) + ' Receiver of No.' + str(lineid) + ' Line (' + \
          linenm + ')')
plt.xlim([np.min(seismot[recid,:]),np.max(seismot[recid,:])])
if flag_figsave:
    subprocess.call('mkdir -p {}'.format(figpath),shell=True)
    figfullnm_rec=figpath + '/' + fignm_rec
    plt.savefig(figfullnm_rec)

# plot receiver line
scl=np.max(np.abs(seismodata))
# increment of show for Y-axis ticklabel
ytickincre=9
#plt.figure()
plt.figure(dpi=figdpi,figsize=(figsize_line[0],figsize_line[1]))
for irec in range(nrec):
    plt.plot(seismot[irec,:],seismodata[irec,:]+irec*(2*scl),'b')
plt.xlabel('Time (s)')
plt.ylabel('Receiver ID')
plt.title(varnm + ' at No.' + str(lineid) + ' Receiver Line (' + linenm + ')')
plt.xlim([np.min(seismot[-1,:]),np.max(seismot[-1,:])])
plt.ylim([-(1+0.5)*(2*scl),(nrec+0.5)*(2*scl)])
plt.yticks(np.arange(0,nrec,ytickincre)*(2*scl),np.arange(1,nrec+1,ytickincre))
if flag_figsave:
    subprocess.call('mkdir -p {}'.format(figpath),shell=True)
    figfullnm_line=figpath + '/' + fignm_line
    plt.savefig(figfullnm_line)

if flag_show:
    plt.show()



