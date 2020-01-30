import numpy as np
import astropy.units as u
from ReadFile import read



def ParticleInfo(filename, par_type, par_num):

    time, total, data = read(filename)

    #x = data['x'][par_num]
    #y = data['y'][par_num]
    #z = data['z'][par_num]
    #vx = data['vx'][par_num]
    #vy = data['vy'][par_num]
    #vz = data['vz'][par_num]
    #m = data['m'][par_num]

    index = np.where(data['type'] == par_type)
     
    xnew = data['x'][index][par_num]
    ynew = data['y'][index][par_num]
    znew = data['z'][index][par_num]
    vxnew = data['vx'][index][par_num]
    vynew = data['vy'][index][par_num]
    vznew = data['vz'][index][par_num]
    mnew = data['m'][index][par_num]
    

    R = np.sqrt(xnew**2 + ynew**2 + znew**2)*u.kpc
    V = np.sqrt(vxnew**2 + vynew**2 + vznew**2)*u.km/u.s
    M = mnew*u.solMass

    return R, V, M

r, v, m = ParticleInfo("MW_000.txt", 2, 99)

#print(r)
print(r.to_value(u.lyr)*u.lyr)
print(v)
print(m)


    
