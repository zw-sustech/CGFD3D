'''
Draw the snapshot on a slice by using pcolor style
Author:      Yuanhang Huo
Email:       yhhuo@mail.ustc.edu.cn
Affiliation: University of Science and Technology of China
Date:        2021.06.25
'''

import numpy as np
import matplotlib.pyplot as plt
import imageio
import subprocess
import json
import glob
import sys
sys.path.append(".")
from netCDF4 import Dataset

# ---------------------------- parameters input --------------------------- #
# file and path name
parfnm='./project/test.json'
slice_dir='./project/output'

# which slice to plot
slicedir='x'
sliceid=19

# which variable and time to plot
varnm='Vz'
ns=1
ne=500
nt=50

# figure control parameters
# 1
flag_show    = 1
taut         = 0.5
# 2
flag_imgsave = 0
flag_gifsave = 1
fignm        = 'fd3dsnap_slice.png'
figsize      = [4,4]
figdpi       = 300
# 3
flag_km      = 1
clbtype      = 'seismic'
clbrange     = [None,None]
# ------------------------------------------------------------------------- #



# read parameter file
par=json.load(open(parfnm,'r'))
ni=np.array(par['number_of_total_grid_points_x'])
nj=np.array(par['number_of_total_grid_points_y'])
nk=np.array(par['number_of_total_grid_points_z'])
nproi=np.array(par['number_of_mpiprocs_x'])
nproj=np.array(par['number_of_mpiprocs_y'])

# slice snapshot plot
plt.figure()

# load data
for nlayer in range(ns,ne+nt,nt):

    # ------------------------ slice x ------------------------ #
    if slicedir == 'x':

        for jp in range(nproj):

            # snapshot data
            slicenm=glob.glob(slice_dir + '/' + 'slicex_i' + str(sliceid) + '_px*_py' + \
                              str(jp) + '.nc')[0]
            slicenc=Dataset(slicenm)
            pnj=slicenc.dimensions['j'].size
            pnk=slicenc.dimensions['k'].size
            if jp == 0:
                V=np.array(slicenc.variables[varnm][nlayer-1,:,:],dtype=np.float64)
            else:
                V0=np.array(slicenc.variables[varnm][nlayer-1,:,:],dtype=np.float64)
                V=np.column_stack((V,V0))
            t=np.array(slicenc.variables['time'][nlayer-1])

            # coordinate data
            slicenm=slicenm.replace(slice_dir + "/","")
            ip=int(slicenm[slicenm.find('px')+2 : slicenm.find('_py')])
            coordnm=slice_dir + '/' + 'coord' + '_px' + str(ip) + '_py' + str(jp) + '.nc'
            coordnc=Dataset(coordnm)
            idwithghost=slicenc.i_index_with_ghosts_in_this_thread
            ghostp=int((coordnc.dimensions['k'].size-slicenc.dimensions['k'].size)/2)
            if jp == 0:
                X=coordnc.variables['x'][ghostp:ghostp+pnk,ghostp:ghostp+pnj,idwithghost]
                Y=coordnc.variables['y'][ghostp:ghostp+pnk,ghostp:ghostp+pnj,idwithghost]
                Z=coordnc.variables['z'][ghostp:ghostp+pnk,ghostp:ghostp+pnj,idwithghost]
            else:
                X0=coordnc.variables['x'][ghostp:ghostp+pnk,ghostp:ghostp+pnj,idwithghost]
                Y0=coordnc.variables['y'][ghostp:ghostp+pnk,ghostp:ghostp+pnj,idwithghost]
                Z0=coordnc.variables['z'][ghostp:ghostp+pnk,ghostp:ghostp+pnj,idwithghost]
                X=np.column_stack((X,X0))
                Y=np.column_stack((Y,Y0))
                Z=np.column_stack((Z,Z0))

        print(V[20,15])
        # units
        str_unit='m'
        if flag_km:
            X=X/1e3
            Y=Y/1e3
            Z=Z/1e3
            str_unit='km'

        # show
        print("Drawing " + str(nlayer) + "th time step (t = " + "{:.4f}".format(t) + "s")
        plt.clf()
        plt.cla()
        plt.pcolor(Y,Z,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])

        plt.xlabel('Y ' + '(' + str_unit + ')')
        plt.ylabel('Z ' + '(' + str_unit + ')')
        plt.title("Snapshot of " + varnm + " at" + "{:9.4f}".format(t) + "s")
        plt.colorbar()
        plt.axis('image')
        if flag_show:
            plt.pause(taut)


if flag_show:
    plt.show()
                









































