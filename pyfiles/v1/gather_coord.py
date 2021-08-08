import numpy as np
import os
from netCDF4 import Dataset

def gather_coord(coordinfo,coorddirvar,coord_dir):

    '''
    gather coordinates data
    Author: Yuanhang Huo
    Email:  yhhuo@mail.ustc.edu.cn
    Date:   2021.06.23
    '''

    if os.path.exists(coord_dir):

        coordprefix='coord'
        nthd=len(coordinfo)

        for n in range(nthd):

            n_i=coordinfo[n]['thisid'][0]
            n_j=coordinfo[n]['thisid'][1]

            i1=coordinfo[n]['indxs'][0]
            j1=coordinfo[n]['indxs'][1]
            k1=coordinfo[n]['indxs'][2]
            
            i2=coordinfo[n]['indxe'][0]
            j2=coordinfo[n]['indxe'][1]
            k2=coordinfo[n]['indxe'][2]

            subs=coordinfo[n]['wsubs']
            subc=coordinfo[n]['wsubc']
            subt=coordinfo[n]['wsubt']

            fnm_coord=coord_dir + "/" + coordprefix  + "_px" + str(n_i) + \
                     "_py" + str(n_j) + ".nc"

            if os.path.exists(fnm_coord):

                coordnc=Dataset(fnm_coord)

                subs=subs[::-1]
                subc=subc[::-1]
                subt=subt[::-1]
                
                xblk=np.array(coordnc.variables['x'][\
                        subs[0]-1:subs[0]+subc[0]*subt[0]-1:subt[0],\
                        subs[1]-1:subs[1]+subc[1]*subt[1]-1:subt[1],\
                        subs[2]-1:subs[2]+subc[2]*subt[2]-1:subt[2]],\
                        dtype=np.float64)

                yblk=np.array(coordnc.variables['y'][\
                        subs[0]-1:subs[0]+subc[0]*subt[0]-1:subt[0],\
                        subs[1]-1:subs[1]+subc[1]*subt[1]-1:subt[1],\
                        subs[2]-1:subs[2]+subc[2]*subt[2]-1:subt[2]],\
                        dtype=np.float64)

                zblk=np.array(coordnc.variables['z'][\
                        subs[0]-1:subs[0]+subc[0]*subt[0]-1:subt[0],\
                        subs[1]-1:subs[1]+subc[1]*subt[1]-1:subt[1],\
                        subs[2]-1:subs[2]+subc[2]*subt[2]-1:subt[2]],\
                        dtype=np.float64)

                if n == 0:
                    xold=xblk
                    yold=yblk
                    zold=zblk

                xnew=np.empty((max(k2,xold.shape[0]),max(j2,xold.shape[1]),max(i2,xold.shape[2])))
                ynew=np.empty((max(k2,yold.shape[0]),max(j2,yold.shape[1]),max(i2,yold.shape[2])))
                znew=np.empty((max(k2,zold.shape[0]),max(j2,zold.shape[1]),max(i2,zold.shape[2])))

                xnew[:xold.shape[0],:xold.shape[1],:xold.shape[2]]=xold
                ynew[:yold.shape[0],:yold.shape[1],:yold.shape[2]]=yold
                znew[:zold.shape[0],:zold.shape[1],:zold.shape[2]]=zold

                xnew[k1-1:k2,j1-1:j2,i1-1:i2]=xblk
                ynew[k1-1:k2,j1-1:j2,i1-1:i2]=yblk
                znew[k1-1:k2,j1-1:j2,i1-1:i2]=zblk

                xold=xnew
                yold=ynew
                zold=znew

            else:
                print("Error in 'gather_coord': File " + "'" + fnm_coord + "'" + \
                      " does NOT exist.")

        xnew=xnew.transpose(2,1,0)
        ynew=ynew.transpose(2,1,0)
        znew=znew.transpose(2,1,0)

        return(xnew,ynew,znew)

    else:
        print("Error in 'gather_coord': Directory " + "'" + coord_dir + "'" + \
              " does NOT exist.")

