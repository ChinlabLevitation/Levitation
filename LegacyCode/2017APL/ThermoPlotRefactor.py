import numpy as np
import math
import matplotlib.pyplot as plt
import os
import thermophoreticForce as thermo
import HeightSimulation as theory

## Plate spacing 6-23: 10.36 mm
## 8-5: 10
os.chdir("C:/Users/Admin/Desktop/Levitating Ice Particles/code/SimulationData/6-23-16-03")

particleRadius = 13.3*10**-6
#particleWeight = 4*math.pi/3*particleRadius**3*998*9.81
#pressure = 6


temperature = []
expHeight = []
theoryHeight = []

for file in os.listdir(os.getcwd()):
    if file.startswith("sim") and file.endswith(".csv"):
        temp = theory.getTheoryHeightArray(file,particleRadius)
        temperature.append(temp[0])
        expHeight.append(1.0-temp[1])
        theoryHeight.append(temp[2])
temperature = np.array(temperature)   
expHeight = np.array(expHeight) 
theoryHeight = np.array(theoryHeight)      

path = "C:/Users/Admin/Desktop/6-23-16-03.txt"
np.savetxt(path,np.transpose([temperature,expHeight]),fmt='%.4f')

plt.scatter(temperature,expHeight*10.0,label='11.3 micron particle',color='black')
plt.xlabel(r'Top Plate Temperature  $\mathit{T}$  [kelvin]',fontsize = 16)
plt.ylabel(r'Particle Height  $\mathit{h}$  [mm]', fontsize = 16)
#plt.scatter(temperature,theoryHeight,c='red',label='Theory')
plt.legend(loc='lower left')
plt.axis([103, 170, 4.7, 7.8])
plt.show()