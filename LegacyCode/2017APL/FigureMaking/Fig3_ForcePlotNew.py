import numpy as np
import math
import matplotlib.pyplot as plt
import os
import thermophoreticForce as thermo
import HeightSimulation as theory
import matplotlib as mpl

mpl.rc('font',family='Times New Roman')
mpl.rcParams['xtick.labelsize'] = 21
mpl.rcParams['ytick.labelsize'] = 21 

plt.subplot(111)
#Set directory where the file is 6-23-16-03/
os.chdir("C:/Users/Real/Desktop/Levitating Ice Particles/code/")
#Load data and initialize values
data = np.loadtxt("C:/Users/Real/Desktop/Levitating Ice Particles/code/SimulationData/tempdata.txt",delimiter=" ")
normalizedData = data[::-1]
normalizedData = theory.cleanData(normalizedData)

heights = np.linspace(0.0,1.0,600)
particleRadius = 13*10**-6
particleWeight = 4*math.pi/3*particleRadius**3*998*9.81
weights = np.full((len(heights)),particleWeight)

plateSpacing = .010
pressure = 6
#Compute thermophoretic force
thermoForce = thermo.calculateForce(normalizedData,particleRadius,plateSpacing,heights,pressure)

#Choose which plots you want
#plt.plot(heights,weights*10**12,c='k',linestyle='--')
#plt.plot(heights*10,thermoForce*10**12,c='k')
#plt.plot(heights*10,-thermo.calculateGradient(normalizedData,heights,plateSpacing)/1000.0,c='k')
#plt.plot(heights*10,thermo.temperature(normalizedData,heights),c='k')

#Force
#plt.xlabel('Height $\mathit{h}$ [mm]',fontsize=18)
#plt.ylabel('$\mathit{F_{th}}$ [pN]',fontsize=18)
#plt.yticks([0,40,80,120,160])

#Temperature Gradient
#plt.xlabel('Height $\mathit{h}$ [mm]',fontsize=18)
#plt.ylabel('Temperature Gradient  $\mathit{T\'}$  [K/mm]',fontsize=18)
#plt.yticks([10,20,30,40])

#Temperature
#plt.xlabel('Height $\mathit{h}$ [mm]',fontsize=18)
#plt.ylabel('Temperature $\mathit{T}$ [K]',fontsize=18)
#plt.yticks([77,177,277])

plt.subplot(3,1,1)
plt.plot(heights*10,thermo.temperature(normalizedData,heights),c='k',linewidth=2)
#plt.ylabel('Temperature $\mathit{T}$ [K]',fontsize=18)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.yticks([77,177,277])
#plt.ylabel('Temperature $\mathit{T}$ [K]',fontsize=18)
plt.text(5.5,220,'Temperature $\mathit{T}$ [K]',fontsize=21)
plt.text(1,100,r'$-\frac{dT}{dz}$',fontsize=28)
plt.text(3,100,'[K/mm]',fontsize=21)
# plt.text(2,80,'Temperature Gradient T\' [K/mm]',fontsize=18)

ax = plt.subplot(3,1,2)
ax = plt.plot(heights*10,-thermo.calculateGradient(normalizedData,heights,plateSpacing)/1000.0,c='k',linewidth=2)
#ax = plt.ylabel('Temperature Gradient  $\mathit{T\'}$  [K/mm]',fontsize=18)
ax = plt.yticks([15,25,35])
ax = plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
#plt.ylabel(r'$\nabla T$  [K/mm]',fontsize=19)

plt.subplot(3,1,3)
plt.plot(heights*10,thermoForce*10**12,c='k',linewidth=2)
plt.plot(heights*10,weights*10**12,c='k',linestyle='--')
plt.yticks([140,90,40])
plt.ylabel('$\mathit{F_{th}}$ [pN]',fontsize=21)
plt.xlabel('Height $\mathit{z}$ [mm]',fontsize=21)

plt.tight_layout(h_pad=0)
plt.subplots_adjust(bottom=.12)

plt.show()