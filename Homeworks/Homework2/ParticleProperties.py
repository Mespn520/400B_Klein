import numpy as np
import astropy.units as u
from ReadFile import read
#A python script that will show values for distance, velocity, and mass for MW_000.txt values

def ParticleInfo(filename, par_type, par_num): #this defines the function ParticleInfor

    time, total, data = read(filename) #This will read the file specifically the values for time, total, and data

    index = np.where(data['type'] == par_type) #Index for all particles with given property of par_type
     
    xnew = data['x'][index][par_num] #the new x-value
    ynew = data['y'][index][par_num] # the new y-value
    znew = data['z'][index][par_num] # the new z-value
    vxnew = data['vx'][index][par_num] # the new vx-value
    vynew = data['vy'][index][par_num] # the new vy-value
    vznew = data['vz'][index][par_num] # the new vz-value
    mnew = data['m'][index][par_num] # the new m-value
    

    R = np.sqrt(xnew**2 + ynew**2 + znew**2)*u.kpc # Radius or magnitude of the distance in kpc
    V = np.sqrt(vxnew**2 + vynew**2 + vznew**2)*u.km/u.s # Velocity or magnitude of the velocity in km/s 
    M = mnew*u.solMass # Mass in units of solar mass

    return R, V, M

r, v, m = ParticleInfo("MW_000.txt", 2, 99) # Function that places "values" on the variables r,v, and m that correspond to what is wanted in the homework

#print(r) # This is available if you want the distance in kpc
print(r.to_value(u.lyr)*u.lyr) # This will convert to ly
print(v) # this is the velocity in km/s
print(m*1e10) # this is the mass in solar masses


    
