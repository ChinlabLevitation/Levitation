import os
os.chdir("C:/Users/Real/Desktop/Levitating Ice Particles/Code/")
import numpy as np
import matplotlib.pyplot as plt
import ftSimulation as ftSim
import matplotlib as mpl

mpl.rc('font',family='Times New Roman')
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20 

plt.subplot(111)

os.chdir("C:/Users/Real/Desktop/Levitating Ice Particles/Code/")
path = "C:/Users/Real/Desktop/PaperData/Data/"

t1 = np.loadtxt(path+"8-5-16_height.txt",usecols=(0,),skiprows=1)
h1 = np.loadtxt(path+"8-5-16_height.txt",usecols=(1,),skiprows=1)

t2 = np.loadtxt(path+"6-23-16-03_height.txt",usecols=(0,),skiprows=1)
h2 = np.loadtxt(path+"6-23-16-03_height.txt",usecols=(1,),skiprows=1)

plt.scatter(t1,h1*10.0,label='radius = 11.3 micron',c='black',s=16,marker='o',facecolor='none')
plt.scatter(t2,h2*10.36,label='radius = 13.3 micron',c='red',s=16,marker='o')
#plt.legend(loc='upper left')
plt.xlabel(r'Top plate temperature  $\mathit{T}$  [K]',fontsize = 21)
plt.ylabel(r'Height  $\mathit{z}$  [mm]', fontsize = 21)
plt.axis([77,180,0,10.4])

mpl.rc('font',family='Times New Roman')

plt.subplots_adjust(bottom=.12)

plt.show()