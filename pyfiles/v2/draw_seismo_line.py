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
import argparse
import sys
sys.path.append(".")
from obspy import read

# ---------------------------- parameters input --------------------------- #
parin=argparse.ArgumentParser(description='Introduction to the drawing script')
parin.add_argument('--parfnm',type=str,required=True,help='parameter json filename with path')
parin.add_argument('--seismo_dir',type=str,required=True,help='seismogram nc file path')
parin.add_argument('--lineid',type=int,required=True,help='receiver line id to plot (1 is the first line))')
parin.add_argument('--recid',type=int,required=True,help='receiver id of the line to plot (0 is the first receiver))')
parin.add_argument('--varnm',type=str,required=True,help="variable to plot, e.g., 'Vp','Tzz'")
parin.add_argument('--flag_show',type=int,default=1,help='show seismogram or not, default=1')
parin.add_argument('--flag_figsave',type=int,default=1,help='save seismo figure or not, default=1')
parin.add_argument('--figpath',type=str,default='./fig',help='figure path to save, default=./fig')
parin.add_argument('--fignm_rec',type=str,default='seismo_rec.png',help='figure name of receiver seismo to save, default=seismo_rec.png')
parin.add_argument('--figsize_rec',type=str,default='[8,4]',help='figure size of receiver seismo to save, default=[8,4]')
parin.add_argument('--fignm_line',type=str,default='seismo_line.png',help='figure name of line seismo to save, default=seismo_line.png')
parin.add_argument('--figsize_line',type=str,default='[5,10]',help='figure size of line seismo to save, default=[5,10]')
parin.add_argument('--figdpi',type=int,default=300,help='figure resolution to save, default=300')

# get all input parameters
par=parin.parse_args()

# parameter json filename with path
parfnm=par.parfnm
# media nc file path
seismo_dir=par.seismo_dir

# which line to plot (1 is the first line)
lineid=par.lineid
# which receiver in the line to plot (0 is the first receiver)
recid=par.recid

# variable name to plot
varnm=par.varnm

# show figure or not
flag_show=par.flag_show
# save figure or not
flag_figsave=par.flag_figsave
# figure path to save
figpath=par.figpath
# figure name of receiver seismo to save
fignm_rec=par.fignm_rec
# figure size of receiver seismo to save
figsize_recstr=par.figsize_rec.split(',')
figsize_rec=[int(figsize_recstr[0][1:]),int(figsize_recstr[1][:-1])]
# figure name of line seismo to save
fignm_line=par.fignm_line
# figure size of line seismo to save
figsize_linestr=par.figsize_line.split(',')
figsize_line=[int(figsize_linestr[0][1:]),int(figsize_linestr[1][:-1])]
# figure resolution to save
figdpi=par.figdpi

## print input parameters for QC
#print(parfnm,type(parfnm))
#print(seismo_dir,type(seismo_dir))
#print(lineid,type(lineid))
#print(recid,type(recid))
#print(varnm,type(varnm))
#print(flag_show,type(flag_show))
#print(flag_figsave,type(flag_figsave))
#print(figpath,type(figpath))
#print(fignm_rec,type(fignm_rec))
#print(figsize_rec,type(figsize_rec))
#print(fignm_line,type(fignm_line))
#print(figsize_line,type(figsize_line))
#print(figdpi,type(figdpi))
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



