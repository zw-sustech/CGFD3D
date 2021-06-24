import numpy as np
import os
import json
import glob
from netCDF4 import Dataset

def locate_media(parfnm,startvar,gsubs,countvar,gsubc,stridevar,gsubt,mediadirvar,media_dir):

    '''
    locate media index in mpi threads
    Author: Yuanhang Huo
    Email:  yhhuo@mail.ustc.edu.cn
    Date:   2021.06.23
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
        media_subs=[1,1,1]
        media_subc=[-1,-1,-1]
        media_subt=[1,1,1]
        media_subs=np.array(media_subs)
        media_subc=np.array(media_subc)
        media_subt=np.array(media_subt)
        ngijk=np.array([par['number_of_total_grid_points_x'],\
                        par['number_of_total_grid_points_y'],\
                        par['number_of_total_grid_points_z']])

        # reset count=-1 to total number
        indx=np.argwhere(media_subc==-1).flatten()
        media_subc[indx]=np.array(np.floor((ngijk[indx]-media_subs[indx])/media_subt[indx])+1,dtype=np.int64)
        media_sube=media_subs+(media_subc-1)*media_subt

        indx=np.argwhere(gsubc==-1).flatten()
        gsubc[indx]=np.array(np.floor((media_subc[indx]-gsubs[indx])/gsubt[indx])+1,dtype=np.int64)
        gsube=gsubs+(gsubc-1)*gsubt

        # search the nc file headers to locate the threads/processors
        mediaprefix='media'
        px=[]
        py=[]
        for mediafullnm in glob.glob(media_dir + "/" + mediaprefix + "*.nc"):
            medianc=Dataset(mediafullnm)
            medianm=mediafullnm.replace(media_dir + "/","")
            xyzs=medianc.global_index_of_first_physical_points
            xs=xyzs[0]
            ys=xyzs[1]
            xyzc=medianc.count_of_physical_points
            xc=xyzc[0]
            yc=xyzc[1]
            xarray=np.arange(xs,xs+xc)
            yarray=np.arange(ys,ys+yc)
            if len(np.intersect1d(np.argwhere(xarray>=gsubs[0]-1),np.argwhere(xarray<=gsube[0]-1))) != 0 \
               and \
               len(np.intersect1d(np.argwhere(yarray>=gsubs[1]-1),np.argwhere(yarray<=gsube[1]-1))) != 0:
                pxtmp=int(medianm[medianm.find('px')+2:medianm.find('_py')])
                px.append(pxtmp)
                pytmp=int(medianm[medianm.find('py')+2:medianm.find('.nc')])
                py.append(pytmp)

        px=np.array(px)
        py=np.array(py)

        # retrieve the coordinate information
        mediainfo=[]
        for ip in np.arange(len(px)):
            subinfo={}
            mediafullnm=media_dir + "/" + mediaprefix + "_px" + str(px[ip]) + "_py" + str(py[ip]) + ".nc"
            medianc=Dataset(mediafullnm)

            xyzs=medianc.global_index_of_first_physical_points
            xs=xyzs[0]
            ys=xyzs[1]
            zs=xyzs[2]
            xyzc=medianc.count_of_physical_points
            xc=xyzc[0]
            yc=xyzc[1]
            zc=xyzc[2]
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

            subinfo['wsubs']=medianc.local_index_of_first_physical_points + \
                             (subinfo['subs']-1)*media_subt+1
            subinfo['wsubc']=subinfo['indxc']
            subinfo['wsubt']=media_subt*gsubt

            subinfo['fnmprefix']=mediaprefix

            subinfo['ttriple']=np.array([gtstart,gtcount,gtstride])
            
            mediainfo.append(subinfo)
        
        return(mediainfo)
    
    else:
        print("Error in 'locate_media': " + "'" + parfnm + "'" + " does NOT exist. Please check it.")
