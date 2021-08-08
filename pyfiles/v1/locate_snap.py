import numpy as np
import os
import json
import glob
from netCDF4 import Dataset

def locate_snap(parfnm,id,startvar,gsubs,countvar,gsubc,stridevar,gsubt,snapdirvar,snap_dir):

    '''
    locate snapshot index in mpi threads
    Author: Yuanhang Huo
    Email:  yhhuo@mail.ustc.edu.cn
    Date:   2021.06.22
    '''

    gtstart  =  1
    gtcount  = -1
    gtstride =  1

    if len(gsubs) == 4:
        gtstart=gsubs[3]
        gsubs=gsubs[0:3]
    if len(gsubc) == 4:
        gtcount=gsubc[3]
        gsubc=gsubc[0:3]
    if len(gsubt) == 4:
        gtstride=gsubt[3]
        gsubt=gsubt[0:3]

    gsubs=np.array(gsubs)
    gsubc=np.array(gsubc)
    gsubt=np.array(gsubt)

    if os.path.exists(parfnm):

        # read parameters file
        par=json.load(open(parfnm,'r'))
        snap_subs=np.array(par['snapshot'][id-1]['grid_index_start'])+1
        snap_subc=np.array(par['snapshot'][id-1]['grid_index_count'])
        snap_subt=np.array(par['snapshot'][id-1]['grid_index_incre'])
        snap_tinv=np.array(par['snapshot'][id-1]['time_index_incre'])
        ngijk=np.array([par['number_of_total_grid_points_x'],\
                        par['number_of_total_grid_points_y'],\
                        par['number_of_total_grid_points_z']])

        if 'snap_subs' in locals():

            # reset count=-1 to total number
            indx=np.argwhere(snap_subc==-1).flatten()
            snap_subc[indx]=np.array(np.floor((ngijk[indx]-snap_subs[indx])/snap_subt[indx])+1,dtype=np.int64)
            snap_sube=snap_subs+(snap_subc-1)*snap_subt

            indx=np.argwhere(gsubc==-1).flatten()
            gsubc[indx]=np.array(np.floor((snap_subc[indx]-gsubs[indx])/gsubt[indx])+1,dtype=np.int64)
            gsube=gsubs+(gsubc-1)*gsubt

            # search the nc file headers to locate the threads/processors
            snapprefix=par['snapshot'][id-1]['name']
            px=[]
            py=[]
            for snapfullnm in glob.glob(snap_dir + "/" + snapprefix + "*.nc"):
                snapnc=Dataset(snapfullnm)
                snapnm=snapfullnm.replace(snap_dir + "/","")
                xyzs=snapnc.first_index_to_snapshot_output
                xs=xyzs[0]
                ys=xyzs[1]
                xc=snapnc.dimensions['i'].size
                yc=snapnc.dimensions['j'].size
                xarray=np.arange(xs,xs+xc)
                yarray=np.arange(ys,ys+yc)
                if len(np.intersect1d(np.argwhere(xarray>=gsubs[0]-1),np.argwhere(xarray<=gsube[0]-1))) != 0 \
                   and \
                   len(np.intersect1d(np.argwhere(yarray>=gsubs[1]-1),np.argwhere(yarray<=gsube[1]-1))) != 0:
                    pxtmp=int(snapnm[snapnm.find('px')+2:snapnm.find('_py')])
                    px.append(pxtmp)
                    pytmp=int(snapnm[snapnm.find('py')+2:snapnm.find('.nc')])
                    py.append(pytmp)

            px=np.array(px)
            py=np.array(py)

            # retrieve the snapshot information
            snapinfo=[]
            for ip in np.arange(len(px)):
                subinfo={}
                snapfullnm=snap_dir + "/" + snapprefix + "_px" + str(px[ip]) + "_py" + str(py[ip]) + ".nc"
                snapnc=Dataset(snapfullnm)

                xyzs=snapnc.first_index_to_snapshot_output
                xs=xyzs[0]
                ys=xyzs[1]
                zs=xyzs[2]
                xc=snapnc.dimensions['i'].size
                yc=snapnc.dimensions['j'].size
                zc=snapnc.dimensions['k'].size
                xe=xs+xc-1
                ye=ys+yc-1
                ze=zs+zc-1

                gxarray=np.arange(gsubs[0],gsube[0]+gsubt[0],gsubt[0])
                gxarray=gxarray-1
                gyarray=np.arange(gsubs[1],gsube[1]+gsubt[1],gsubt[1])
                gyarray=gyarray-1
                gzarray=np.arange(gsubs[2],gsube[2]+gsubt[2],gsubt[2])
                gzarray=gzarray-1

                subinfo['thisid']=np.array([px[ip],py[ip]])
                i=np.intersect1d(np.argwhere(gxarray>=xs),np.argwhere(gxarray<=xe))+1
                j=np.intersect1d(np.argwhere(gyarray>=ys),np.argwhere(gyarray<=ye))+1
                k=np.intersect1d(np.argwhere(gzarray>=zs),np.argwhere(gzarray<=ze))+1
                subinfo['indxs']=np.array([i[0],j[0],k[0]])
                subinfo['indxe']=np.array([i[-1],j[-1],k[-1]])
                subinfo['indxc']=subinfo['indxe']-subinfo['indxs']+1

                subinfo['subs']=np.array([gxarray[i[0]-1]-xs+1,\
                                          gyarray[j[0]-1]-ys+1,\
                                          gzarray[k[0]-1]-zs+1])
                subinfo['subc']=subinfo['indxc']
                subinfo['subt']=gsubt

                subinfo['wsubs']=snapnc.first_index_in_this_thread_with_ghosts + \
                                 (subinfo['subs']-1)*snap_subt+1
                subinfo['wsubc']=subinfo['indxc']
                subinfo['wsubt']=snap_subt*gsubt

                subinfo['tinv']=snap_tinv

                subinfo['fnmprefix']=snapprefix

                subinfo['ttriple']=np.array([gtstart,gtcount,gtstride])
                
                snapinfo.append(subinfo)
            
            return(snapinfo)
    
        else:
            print("Error in 'locate_snap': " + "id = " + str(id) + " does NOT exist.")

    else:
        print("Error in 'locate_snap': " + "'" + parfnm + "'" + " does NOT exist. Please check it.")
