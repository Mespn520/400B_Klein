import numpy as np
import astropy.units as u
#A python script that will read and use the values defined in the MW_000.txt file to evaluate different things.

filename = "MW_000.txt" #this defines the filename

def read(filename): #This defines the function read that takes the name of the file as input
    
    file = open(filename, 'r') #This will open the file
    line1 = file.readline() #This is line1
    label, value = line1.split() #This is how it is labelled and valued
    time = float(value)*u.Myr #This will store the time in units of Myr
    
    line2 = file.readline() #This is line2
    label, value = line2.split() #This is how it is labelled and valued
    total = float(value) #This will store the total particles
    file.close() #closes the file
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3) #This allows for the use of the column header information that is presented in the txt file starting with a #-symbol
    #print(data['type'][n]) #n just defines the row number 
    #print(data['x'][n])
    #print(data['y'][n])
    #print(data['z'][n])
    #print(data['vx'][n])
    #print(data['vy'][n])
    #print(data['vz'][n])
    #print(data['m'][n])
    
    return time,total,data #This returns the time, total particles and data

if __name__ == "__main__":
    time, total, data = read(filename) #This just reads the file

    print("The time is",time) #This prints the time
    print("The total particles are",total) #This prints the total particles 
    print(data) #This prints the data array
    
    

