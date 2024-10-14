import numpy as np
import math
import matplotlib.pyplot as plt
import os
import thermophoreticForce as thermo
import HeightSimulation as theory



os.chdir("C:/Users/Admin/Desktop/Levitating Ice Particles/code/SimulationData/6-23-16-03")

particleRadius = 12*10**-6
#particleWeight = 4*math.pi/3*particleRadius**3*998*9.81
#pressure = 6

kn = np.array([])
temperature = np.array([])
expHeight = []
theoryHeight = []

for file in os.listdir(os.getcwd()):
    if file.startswith("sim") and file.endswith(".csv"):
            print(file)
            data =theory.getData(file)
            data = theory.cleanData(data)
            constants = theory.getConstants(file)
    
            plateSpacing = constants[2]/float(1000.0)
            pressure = constants[3]
            experimentalHeight = 1.0-constants[4] 
            topPlateTemp = data[len(data)-1]
            kn = np.append(kn,thermo.computeKnudsenNumber(data,experimentalHeight,pressure,particleRadius))
            temperature = np.append(temperature,topPlateTemp)

#temperature = np.array(temperature)   
#expHeight = np.array(expHeight) 
#theoryHeight = np.array(theoryHeight)      

plt.scatter(temperature,kn,label='Data')
#plt.scatter(temperature,theoryHeight,c='red',label='Theory')
#plt.legend()
#plt.axis([80, 170, -.1, 1.1])
plt.show()