import numpy as np
import astropy.units as u
from ReadFile import read
import matplotlib.pyplot as plt
import matplotlib
from CenterOfMass import CenterOfMass
from astropy.constants import G

filename1="./MW_000.txt"
filename11="./M31_000.txt"
filename111="./M33_000.txt"

data1 = np.genfromtxt(filename1,dtype=None,names=True,skip_header=4)
data11 = np.genfromtxt(filename11,dtype=None,names=True,skip_header=4)
data111 = np.genfromtxt(filename111,dtype=None,names=True,skip_header=4)

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

plt.plot(data1, color='blue', linewidth=5, label='MW')
plt.plot(data11, color='red', linewidth=5, label='M31')
plt.plot(data111, color='green', linewidth=5, label='M33')

plt.xlabel('Galaxy', fontsize=22)
plt.ylabel('M$_dot$', fontsize=22)

plt.xlim(0, 30)
plt.ylim(1e6, 1e12)

label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size

legend = ax.legend(loc='upper left', fontsize='x-large')

plt.figtext(0.6, 0.15, 'Galaxies', fontsize=22)
