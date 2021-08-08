import numpy as np
import os
from netCDF4 import Dataset

def gather_metric(metricinfo,varnm,metricdirvar,metric_dir):

    '''
    gather metric data
    Author: Yuanhang Huo
    Email:  yhhuo@mail.ustc.edu.cn
    Date:   2021.06.24
    '''

    if os.path.exists(metric_dir):

        metricprefix='metric'
        nthd=len(metricinfo)

        for n in range(nthd):

            n_i=metricinfo[n]['thisid'][0]
            n_j=metricinfo[n]['thisid'][1]

            i1=metricinfo[n]['indxs'][0]
            j1=metricinfo[n]['indxs'][1]
            k1=metricinfo[n]['indxs'][2]
            
            i2=metricinfo[n]['indxe'][0]
            j2=metricinfo[n]['indxe'][1]
            k2=metricinfo[n]['indxe'][2]

            subs=metricinfo[n]['wsubs']
            subc=metricinfo[n]['wsubc']
            subt=metricinfo[n]['wsubt']

            fnm_metric=metric_dir + "/" + metricprefix  + "_px" + str(n_i) + \
                     "_py" + str(n_j) + ".nc"

            if os.path.exists(fnm_metric):

                metricnc=Dataset(fnm_metric)

                subs=subs[::-1]
                subc=subc[::-1]
                subt=subt[::-1]
                
                vblk=np.array(metricnc.variables[varnm][\
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
                print("Error in 'gather_metric': File " + "'" + fnm_metric + "'" + \
                      " does NOT exist.")

        vnew=vnew.transpose(2,1,0)

        return(vnew)

    else:
        print("Error in 'gather_metric': Directory " + "'" + metric_dir + "'" + \
              " does NOT exist.")

