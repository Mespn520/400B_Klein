# Homework 5
# Determine the mass distribution of each galaxy at SnapNumber 0 and use this to determine each galaxy's rotation curve.
# Michael Klein
# Worked with Trevor Smith and James Taylor

# import modules
import numpy as np
import astropy.units as u
from ReadFile import read
import matplotlib.pyplot as plt
import matplotlib
from CenterOfMass import CenterOfMass
import astropy.constants as const



class MassProfile:
# Class to define the MassProfile of the galaxies

    def __init__(self, galaxy, snap):
        # add a string of the filenumber to the value "000"
        ilbl = '000' + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy) + ilbl + '.txt'
        
        
        # Reading in the data
        self.time, self.total, self.data = read(self.filename)
        self.type = self.data['type']
        self.m = self.data['m']
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc

        self.gname = galaxy

    def MassEnclosed(self, ptype, r):
        # Inputs:
        # ptype : particle type
        # r : magnitude radius
        # Returns:
        # An array of masses (in Msun) so we can use it to calculate the mass profile conveniently
        # Creating the COM object to then call COM position
        COM = CenterOfMass(self.filename, ptype)

        # Calling COM position
        delta = 0.1
        XCOM, YCOM, ZCOM = COM.COM_P(delta)

        # change reference frame to COM frame

        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM

        index = np.where(self.type == ptype)
        xNew2 = xNew[index]
        yNew2 = yNew[index]
        zNew2 = zNew[index]

        # Calculate magnitude
        R = (xNew2**2 + yNew2**2 + zNew2**2)**0.5
        
        # Create array
        enmass = np.zeros(len(r))
        
        # Loop over the radius array
        for i in range(len(r)):
            index2 = np.where(R < r[i])
            # Adding the masses
            enmass[i] = np.sum(self.m[index2])
            
        return enmass*1e10*u.Msun
    
    def MassEnclosedTotal(self, r):
        # Inputs:
        # r : The array of radii (1D array)
        # Returns:
        # An array of masses (Msun) representing the total enclosed mass (bulge+disk+halo) at each radius of the input array

        # Getting mass for each ptype
        DiskMass = self.MassEnclosed(2,r)
        HaloMass = self.MassEnclosed(1,r)

        # return the sum of these at r
    # If statement to take into account that M33 does not contain a bulge 
        if self.gname == 'M33':
            BulgeMass = np.zeros(len(r))
        else:
            BulgeMass = self.MassEnclosed(3,r)
        return DiskMass + HaloMass + BulgeMass

    def HernquistMass(self,r,a,Mhalo):
        # Inputs:
        # r : The radius (kpc)
        # a : The scale factor (kpc)
        # Mhalo : The halo mass (Msun)
        # Returns:
        # The halo mass in units of Msun
        return (Mhalo*r**2)/((a+r)**2)

    def CircularVelocity(self, ptype, r):
        # Inputs:
        # ptype : The particle type
        # r : The radius (kpc)
        # Returns:
        # An array of circular speeds in units of km/s, rounding to two decimal places
        # Adjusting G to the correct units this program works in
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)

        return np.round(np.sqrt(G*self.MassEnclosed(ptype,r)/r),2)

    def CircularVelocityTotal(self, r):
        # Inputs:
        # ptype : The particle type
        # r : The radius (kpc)
        # Returns:
        # An array of circular speeds in units of km/s, rounding to two decimal places
        # Adjusting G to the correct units this program works in
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)

        return np.round(np.sqrt(G*self.MassEnclosedTotal(r)/r),2)


    def HernquistVCirc(self, r, a, Mhalo):
        # Inputs:
        # r : The radius (kpc)
        # a : The scale factor (kpc)
        # Mhalo : The halo Mass (Msun)
        # Returns:
        # The circular speed in units of km/s, rounded to two decimals

        # Adjust G to correct units
        G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)

        return np.round(np.sqrt(G*self.HernquistMass(r, a, Mhalo)/r),2)
    
# Plotting mass profile for each galaxy
# Setting an array of radii
    
r = np.linspace(0.1, 30)*u.kpc
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(30,20))
MassP_MW = MassProfile('MW',0)
MW_Mass_Halo = MassP_MW.MassEnclosed(1,r)
MW_Mass_Disk = MassP_MW.MassEnclosed(2,r)
MW_Mass_Bulge = MassP_MW.MassEnclosed(3,r)
MassP_M31 = MassProfile('M31',0)
M31_Mass_Halo = MassP_M31.MassEnclosed(1,r)
M31_Mass_Disk = MassP_M31.MassEnclosed(2,r)
M31_Mass_Bulge = MassP_M31.MassEnclosed(3,r)
MassP_M33 = MassProfile('M33',0)
M33_Mass_Halo = MassP_M33.MassEnclosed(1,r)
M33_Mass_Disk = MassP_M33.MassEnclosed(2,r)


ax[0,0].semilogy(r, MW_Mass_Halo, label='Halo Mass', color='black', linewidth=3, linestyle=':')
ax[0,0].semilogy(r, MW_Mass_Disk, label='Disk Mass', color='red' , linewidth=3, linestyle='-')
ax[0,0].semilogy(r, MW_Mass_Bulge, label='Bulge Mass', color='green', linewidth=3,linestyle='--')
ax[0,1].semilogy(r, M31_Mass_Halo, label='Halo Mass', color='black', linewidth=3, linestyle=':')
ax[0,1].semilogy(r, M31_Mass_Disk, label='Disk Mass', color='red' , linewidth=3, linestyle='-')
ax[0,1].semilogy(r, M31_Mass_Bulge, label='Bulge Mass', color='green', linewidth=3,linestyle='--')
ax[0,2].semilogy(r, M33_Mass_Halo, label='Halo Mass', color='black', linewidth=3, linestyle=':')
ax[0,2].semilogy(r, M33_Mass_Disk, label='Disk Mass', color='red' , linewidth=3, linestyle='-')


ax[0,0].set(title='MW',xlabel='Radius (kpc)', ylabel='Mass Enclosed $M_{\odot}$')
ax[0,1].set(title='M31',xlabel='Radius (kpc)', ylabel='Mass Enclosed $M_{\odot}$')
ax[0,2].set(title='M33',xlabel='Radius (kpc)', ylabel='Mass Enclosed $M_{\odot}$')



# Plotting the Hernquist Mass Profiles

indexMW = np.where(MassP_MW.type==1)
MWHalo = np.sum(MassP_MW.m[indexMW])*1e10
indexM31 = np.where(MassP_M31.type==1)
M31Halo = np.sum(MassP_M31.m[indexM31])*1e10
indexM33 = np.where(MassP_M33.type==1)
M33Halo = np.sum(MassP_M33.m[indexM33])*1e10

MW_HernMass = MassP_MW.HernquistMass(r,1*u.kpc,MWHalo)
M31_HernMass = MassP_M31.HernquistMass(r,1.5*u.kpc,M31Halo)
M33_HernMass = MassP_M33.HernquistMass(r,0.5*u.kpc,M33Halo)

ax[0,0].semilogy(r, MW_HernMass, label='MW Hernquist Mass', color='blue' , linewidth=3,linestyle='-.')
ax[0,1].semilogy(r, M31_HernMass, label='M31 Hernquist Mass', color='blue' , linewidth=3,linestyle='-.')
ax[0,2].semilogy(r, M33_HernMass, label='M33 Hernquist Mass', color='blue' , linewidth=3,linestyle='-.')

#ax[0,0].set(title='MW Hernquist Mass',xlabel='Distance (kpc)', ylabel='Mass $M_{\odot}$')
#ax[1,1].set(title='M31 Hernquist Mass',xlabel='Distance (kpc)', ylabel='Mass $M_{\odot}$')
#ax[1,2].set(title='M33 Hernquist Mass',xlabel='Distance (kpc)', ylabel='Mass $M_{\odot}$')

# Plotting Circular Velocity
MW_CircV_Halo = MassP_MW.CircularVelocity(1,r)
MW_CircV_Disk = MassP_MW.CircularVelocity(2,r)
MW_CircV_Bulge = MassP_MW.CircularVelocity(3,r)
M31_CircV_Halo = MassP_M31.CircularVelocity(1,r)
M31_CircV_Disk = MassP_M31.CircularVelocity(2,r)
M31_CircV_Bulge = MassP_M31.CircularVelocity(3,r)
M33_CircV_Halo = MassP_M33.CircularVelocity(1,r)
M33_CircV_Disk = MassP_M33.CircularVelocity(2,r)
MW_CircV_Total = MassP_MW.CircularVelocityTotal(r)
M31_CircV_Total = MassP_M31.CircularVelocityTotal(r)
M33_CircV_Total = MassP_M33.CircularVelocityTotal(r)

ax[1,0].semilogy(r, MW_CircV_Halo, label='MW Halo', color='blue' , linewidth=3,linestyle='-.')
ax[1,0].semilogy(r, MW_CircV_Disk, label='MW Disk', color='green' , linewidth=3,linestyle='--')
ax[1,0].semilogy(r, MW_CircV_Bulge, label='MW Bulge', color='red' , linewidth=3,linestyle=':')
ax[1,0].semilogy(r, MW_CircV_Total, label='MW Total', color='black' , linewidth=3)
ax[1,1].semilogy(r, M31_CircV_Halo, label='M31 Halo', color='blue' , linewidth=3,linestyle='-.')
ax[1,1].semilogy(r, M31_CircV_Disk, label='M31 Disk', color='green' , linewidth=3,linestyle='--')
ax[1,1].semilogy(r, M31_CircV_Bulge, label='M31 Bulge', color='red' , linewidth=3,linestyle=':')
ax[1,1].semilogy(r, M31_CircV_Total, label='M31 Total', color='black' , linewidth=3)
ax[1,2].semilogy(r, M33_CircV_Halo, label='M33 Halo', color='blue' , linewidth=3,linestyle='-.')
ax[1,2].semilogy(r, M33_CircV_Disk, label='M33 Disk', color='green' , linewidth=3,linestyle='--')
ax[1,2].semilogy(r, M33_CircV_Total, label='M33 Total', color='black' , linewidth=3)



ax[0,0].legend()
ax[0,1].legend()
ax[0,2].legend()
ax[1,0].legend()
ax[1,1].legend()
ax[1,2].legend()

#plt.tight_layout()
plt.show()
