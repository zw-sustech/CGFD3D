import numpy as np
import os
from netCDF4 import Dataset

def gather_media(mediainfo,varnm,mediadirvar,media_dir):

    '''
    gather media data
    Author: Yuanhang Huo
    Email:  yhhuo@mail.ustc.edu.cn
    Date:   2021.06.23
    '''

    if os.path.exists(media_dir):

        mediaprefix='media'
        nthd=len(mediainfo)

        for n in range(nthd):

            n_i=mediainfo[n]['thisid'][0]
            n_j=mediainfo[n]['thisid'][1]

            i1=mediainfo[n]['indxs'][0]
            j1=mediainfo[n]['indxs'][1]
            k1=mediainfo[n]['indxs'][2]
            
            i2=mediainfo[n]['indxe'][0]
            j2=mediainfo[n]['indxe'][1]
            k2=mediainfo[n]['indxe'][2]

            subs=mediainfo[n]['wsubs']
            subc=mediainfo[n]['wsubc']
            subt=mediainfo[n]['wsubt']

            fnm_media=media_dir + "/" + mediaprefix  + "_px" + str(n_i) + \
                     "_py" + str(n_j) + ".nc"

            if os.path.exists(fnm_media):

                medianc=Dataset(fnm_media)

                subs=subs[::-1]
                subc=subc[::-1]
                subt=subt[::-1]
                
                vblk=np.array(medianc.variables[varnm][\
                        subs[0]-1:subs[0]+subc[0]*subt[0]-1:subt[0],\
                        subs[1]-1:subs[1]+subc[1]*subt[1]-1:subt[1],\
                        subs[2]-1:subs[2]+subc[2]*subt[2]-1:subt[2]],\
                        dtype=np.float64)

                if n == 0:
                    vold=vblk

                vnew=np.empty((max(k2,vold.shape[0]),max(j2,vold.shape[1]),max(i2,vold.shape[2])))

                vnew[:vold.shape[0],:vold.shape[1],:vold.shape[2]]=vold

                vnew[k1-1:k2,j1-1:j2,i1-1:i2]=vblk

                vold=vnew

            else:
                print("Error in 'gather_media': File " + "'" + fnm_media + "'" + \
                      " does NOT exist.")

        vnew=vnew.transpose(2,1,0)

        return(vnew)

    else:
        print("Error in 'gather_media': Directory " + "'" + media_dir + "'" + \
              " does NOT exist.")

