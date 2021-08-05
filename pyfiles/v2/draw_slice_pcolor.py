'''
Draw the snapshot on a slice by using pcolor style
Author:      Yuanhang Huo
Email:       yhhuo@mail.ustc.edu.cn
Affiliation: University of Science and Technology of China
Date:        2021.06.26
'''

import numpy as np
import matplotlib.pyplot as plt
import imageio
import subprocess
import argparse
import json
import glob
import sys
sys.path.append(".")
from netCDF4 import Dataset

# ---------------------------- parameters input --------------------------- #
parin=argparse.ArgumentParser(description='Introduction to the drawing script')
parin.add_argument('--parfnm',type=str,required=True,help='parameter json filename with path')
parin.add_argument('--slice_dir',type=str,required=True,help='slice snapshot nc file path')
parin.add_argument('--slicedir',type=str,required=True,help='which direction of slice to path')
parin.add_argument('--sliceid',type=int,required=True,help='slice id')
parin.add_argument('--varnm',type=str,required=True,help="variable to plot, e.g., 'Vx','Tzz'")
parin.add_argument('--ns',type=int,required=True,help='starting time step to plot')
parin.add_argument('--ne',type=int,required=True,help='ending time step to plot')
parin.add_argument('--nt',type=int,required=True,help='time stride to plot')
parin.add_argument('--flag_show',type=int,default=1,help='show snapshot or not, default=1')
parin.add_argument('--taut',type=float,default=0.5,help='plotting pause time (in second) between two time steps, default=0.5')
parin.add_argument('--flag_imgsave',type=int,default=1,help='save snapshot figure or not, default=1')
parin.add_argument('--flag_gifsave',type=int,default=1,help='save snapshot gif or not, default=1')
parin.add_argument('--fignm',type=str,default='fd3dsnap_slice.png',help='figure name to save, default=fd3dsnap_slice.png')
parin.add_argument('--figsize',type=str,default='[4,4]',help='figure size to save, default=[4,4]')
parin.add_argument('--figdpi',type=int,default=300,help='figure resolution to save, default=300')
parin.add_argument('--flag_km',type=int,default=1,help='figure axis unit in km or not, default=1')
parin.add_argument('--clbtype',type=str,default='seismic',help='colorbar type, default=seismic')
parin.add_argument('--clbrange',type=str,default='[None,None]',help='colorbar range, default=[None,None]')

# get all input parameters
par=parin.parse_args()

# parameter json filename with path
parfnm=par.parfnm
# slice snapshot nc file path
slice_dir=par.slice_dir

# slice direction to plot
slicedir=par.slicedir
# slice id
sliceid=par.sliceid

# variable name to plot
varnm=par.varnm
# starting time step to plot
ns=par.ns
# ending time step to plot
ne=par.ne
# time stride to plot
nt=par.nt

# show figure or not
flag_show=par.flag_show
# plotting pause time in second
taut=par.taut
# save figure or not
flag_imgsave=par.flag_imgsave
# save gif or not
flag_gifsave=par.flag_gifsave
# figure name to save
fignm=par.fignm
# figure size to save
figsizestr=par.figsize.split(',')
figsize=[int(figsizestr[0][1:]),int(figsizestr[1][:-1])]
# figure resolution to save
figdpi=par.figdpi
# axis unit km or m
flag_km=par.flag_km
# colorbar type
clbtype=par.clbtype
# colorbar range
clbrangestr=par.clbrange.split(',')
if clbrangestr[0][1:] == 'None' and clbrangestr[1][:-1] == 'None':
    clbrange=[None,None]
else:
    clbrange=[float(clbrangestr[0][1:]),float(clbrangestr[1][:-1])]

## print input parameters for QC
#print(parfnm,type(parfnm))
#print(slice_dir,type(slice_dir))
#print(slicedir,type(slicedir))
#print(sliceid,type(sliceid))
#print(varnm,type(varnm))
#print(ns,type(ns))
#print(ne,type(ne))
#print(nt,type(nt))
#print(flag_show,type(flag_show))
#print(taut,type(taut))
#print(flag_imgsave,type(flag_imgsave))
#print(flag_gifsave,type(flag_gifsave))
#print(fignm,type(fignm))
#print(figsize,type(figsize))
#print(figdpi,type(figdpi))
#print(flag_km,type(flag_km))
#print(clbtype,type(clbtype))
#print(clbrange,type(clbrange))
# ------------------------------------------------------------------------- #



# read parameter file
par=json.load(open(parfnm,'r'))
ni=np.array(par['number_of_total_grid_points_x'])
nj=np.array(par['number_of_total_grid_points_y'])
nk=np.array(par['number_of_total_grid_points_z'])
nproi=np.array(par['number_of_mpiprocs_x'])
nproj=np.array(par['number_of_mpiprocs_y'])

# slice snapshot plot
#plt.figure()
plt.figure(dpi=figdpi,figsize=(figsize[0],figsize[1]))

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

    # ------------------------ slice y ------------------------ #
    elif slicedir == 'y':

        for ip in range(nproi):

            # snapshot data
            slicenm=glob.glob(slice_dir + '/' + 'slicey_j' + str(sliceid) + '_px' + \
                              str(ip) + '_py*.nc')[0]
            slicenc=Dataset(slicenm)
            pni=slicenc.dimensions['i'].size
            pnk=slicenc.dimensions['k'].size
            if ip == 0:
                V=np.array(slicenc.variables[varnm][nlayer-1,:,:],dtype=np.float64)
            else:
                V0=np.array(slicenc.variables[varnm][nlayer-1,:,:],dtype=np.float64)
                V=np.column_stack((V,V0))
            t=np.array(slicenc.variables['time'][nlayer-1])

            # coordinate data
            slicenm=slicenm.replace(slice_dir + "/","")
            jp=int(slicenm[slicenm.find('py')+2 : slicenm.find('.nc')])
            coordnm=slice_dir + '/' + 'coord' + '_px' + str(ip) + '_py' + str(jp) + '.nc'
            coordnc=Dataset(coordnm)
            idwithghost=slicenc.j_index_with_ghosts_in_this_thread
            ghostp=int((coordnc.dimensions['i'].size-slicenc.dimensions['i'].size)/2)
            if ip == 0:
                X=coordnc.variables['x'][ghostp:ghostp+pnk,idwithghost,ghostp:ghostp+pni]
                Y=coordnc.variables['y'][ghostp:ghostp+pnk,idwithghost,ghostp:ghostp+pni]
                Z=coordnc.variables['z'][ghostp:ghostp+pnk,idwithghost,ghostp:ghostp+pni]
            else:
                X0=coordnc.variables['x'][ghostp:ghostp+pnk,idwithghost,ghostp:ghostp+pni]
                Y0=coordnc.variables['y'][ghostp:ghostp+pnk,idwithghost,ghostp:ghostp+pni]
                Z0=coordnc.variables['z'][ghostp:ghostp+pnk,idwithghost,ghostp:ghostp+pni]
                X=np.column_stack((X,X0))
                Y=np.column_stack((Y,Y0))
                Z=np.column_stack((Z,Z0))

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
        plt.pcolor(X,Z,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])

        plt.xlabel('X ' + '(' + str_unit + ')')
        plt.ylabel('Z ' + '(' + str_unit + ')')
        plt.title("Snapshot of " + varnm + " at" + "{:9.4f}".format(t) + "s")
        plt.colorbar()
        plt.axis('image')
        if flag_show:
            plt.pause(taut)

    # ------------------------ slice z ------------------------ #
    else:

        for jp in range(nproj):
            for ip in range(nproi):

                # snapshot data
                slicenm=glob.glob(slice_dir + '/' + 'slicez_k' + str(sliceid) + '_px' + \
                                  str(ip) + '_py' + str(jp) + '.nc')[0]
                slicenc=Dataset(slicenm)
                pni=slicenc.dimensions['i'].size
                pnj=slicenc.dimensions['j'].size
                if ip == 0:
                    VV=np.array(slicenc.variables[varnm][nlayer-1,:,:],dtype=np.float64)
                else:
                    VV0=np.array(slicenc.variables[varnm][nlayer-1,:,:],dtype=np.float64)
                    VV=np.column_stack((VV,VV0))
                t=np.array(slicenc.variables['time'][nlayer-1])

                # coordinate data
                coordnm=slice_dir + '/' + 'coord' + '_px' + str(ip) + '_py' + str(jp) + '.nc'
                coordnc=Dataset(coordnm)
                idwithghost=slicenc.k_index_with_ghosts_in_this_thread
                ghostp=int((coordnc.dimensions['j'].size-slicenc.dimensions['j'].size)/2)
                if ip == 0:
                    XX=coordnc.variables['x'][idwithghost,ghostp:ghostp+pnj,ghostp:ghostp+pni]
                    YY=coordnc.variables['y'][idwithghost,ghostp:ghostp+pnj,ghostp:ghostp+pni]
                    ZZ=coordnc.variables['z'][idwithghost,ghostp:ghostp+pnj,ghostp:ghostp+pni]
                else:
                    XX0=coordnc.variables['x'][idwithghost,ghostp:ghostp+pnj,ghostp:ghostp+pni]
                    YY0=coordnc.variables['y'][idwithghost,ghostp:ghostp+pnj,ghostp:ghostp+pni]
                    ZZ0=coordnc.variables['z'][idwithghost,ghostp:ghostp+pnj,ghostp:ghostp+pni]
                    XX=np.column_stack((XX,XX0))
                    YY=np.column_stack((YY,YY0))
                    ZZ=np.column_stack((ZZ,ZZ0))

            if jp == 0:
                V=VV
                X=XX
                Y=YY
                Z=ZZ
            else:
                V=np.row_stack((V,VV))
                X=np.row_stack((X,XX))
                Y=np.row_stack((Y,YY))
                Z=np.row_stack((Z,ZZ))

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
        plt.pcolor(X,Y,V,cmap=clbtype,vmin=clbrange[0],vmax=clbrange[1])

        plt.xlabel('X ' + '(' + str_unit + ')')
        plt.ylabel('Y ' + '(' + str_unit + ')')
        plt.title("Snapshot of " + varnm + " at" + "{:9.4f}".format(t) + "s")
        plt.colorbar()
        plt.axis('image')
        if flag_show:
            plt.pause(taut)

    # save figure and GIF
    if flag_imgsave or flag_gifsave:
        imgnm=fignm[:-(len(fignm.split('.')[-1])+1)]
        imgfmt=fignm.split('.')[-1]
        imgfullnm='{}_timestep_{}.{}'.format(imgnm,nlayer,imgfmt)
        plt.savefig(imgfullnm)

if flag_gifsave:
    frames=[]
    for nlayer in range(ns,ne+nt,nt):
        imgfullnm='{}_timestep_{}.{}'.format(imgnm,nlayer,imgfmt)
        frames.append(imageio.imread(imgfullnm))
    imageio.mimsave('{}.gif'.format(imgnm),frames,'GIF',duration=taut)

if flag_imgsave == 0 and flag_gifsave == 1:
    subprocess.call('rm {}_timestep_*.{}'.format(imgnm,imgfmt),shell=True)


if flag_show:
    plt.show()
                

