import numpy as np
import matplotlib.pyplot as plt
import ftSimulation as ftSim
import math

path = "C:/Users/Admin/Desktop/PaperData/"

t1 = np.loadtxt(path+"8-5-16_height.txt",usecols=(0,),skiprows=1)
h1 = np.loadtxt(path+"8-5-16_height.txt",usecols=(1,),skiprows=1)

t2 = np.loadtxt(path+"6-23-16-03_height.txt",usecols=(0,),skiprows=1)
h2 = np.loadtxt(path+"6-23-16-03_height.txt",usecols=(1,),skiprows=1)

plt.scatter(t1,h1*10.0,label='radius = 11.3 micron',c='black',s=12,marker='o',facecolor='none')
plt.scatter(t2,h2*10.36,label='radius = 13.3 micron',c='black',s=12,marker='o')
#plt.legend(loc='upper left')
plt.xlabel(r'Top plate temperature  $\mathit{T}$  [K]',fontsize = 16)
plt.ylabel(r'Particle height  $\mathit{h}$  [mm]', fontsize = 16)
plt.axis([77,180,0,10.4])
plt.show()