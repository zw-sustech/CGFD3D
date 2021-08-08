import numpy as np
import os
from netCDF4 import Dataset

def gather_snap(snapinfo,nlayer,varnm,snapdirvar,snap_dir):

    '''
    gather snapshot data
    Author: Yuanhang Huo
    Email:  yhhuo@mail.ustc.edu.cn
    Date:   2021.06.23
    '''

    if os.path.exists(snap_dir):

        nthd=len(snapinfo)

        for n in range(nthd):

            n_i=snapinfo[n]['thisid'][0]
            n_j=snapinfo[n]['thisid'][1]

            i1=snapinfo[n]['indxs'][0]
            j1=snapinfo[n]['indxs'][1]
            k1=snapinfo[n]['indxs'][2]
            
            i2=snapinfo[n]['indxe'][0]
            j2=snapinfo[n]['indxe'][1]
            k2=snapinfo[n]['indxe'][2]

            subs=snapinfo[n]['subs']
            subc=snapinfo[n]['subc']
            subt=snapinfo[n]['subt']

            fnm_snap=snap_dir + "/" + snapinfo[n]['fnmprefix'] + "_px" + str(n_i) + \
                     "_py" + str(n_j) + ".nc"

            if os.path.exists(fnm_snap):

                snapnc=Dataset(fnm_snap)
                tdim={}
                tdim['Name']=snapnc.dimensions['time'].name
                tdim['Length']=snapnc.dimensions['time'].size

                if tdim['Length'] == 0 or (nlayer-1)-1 >= tdim['Length']:
                    print("Error in 'gather_snap': " + str(nlayer) + "th layer is beyond \
                           current time dim (" + str(tdim['Length']) + ") in " + fnm_snap)
                else:
                    subs=subs[::-1]
                    subc=subc[::-1]
                    subt=subt[::-1]
                    
                    vblk=np.array(snapnc.variables[varnm][nlayer-1,\
                            subs[0]-1:subs[0]+subc[0]*subt[0]-1:subt[0],\
                            subs[1]-1:subs[1]+subc[1]*subt[1]-1:subt[1],\
                            subs[2]-1:subs[2]+subc[2]*subt[2]-1:subt[2]],\
                            dtype=np.float64)

                    if n == 0:
                        vold=vblk
                        t=snapnc.variables['time'][nlayer-1]

                    vnew=np.empty((max(k2,vold.shape[0]),max(j2,vold.shape[1]),max(i2,vold.shape[2])))

                    vnew[:vold.shape[0],:vold.shape[1],:vold.shape[2]]=vold

                    vnew[k1-1:k2,j1-1:j2,i1-1:i2]=vblk

                    vold=vnew

            else:
                print("Error in 'gather_snap': File " + "'" + fnm_snap + "'" + \
                      " does NOT exist.")

        vnew=vnew.transpose(2,1,0)

        return(vnew,t)

    else:
        print("Error in 'gather_snap': Directory " + "'" + snap_dir + "'" + \
              " does NOT exist.")

