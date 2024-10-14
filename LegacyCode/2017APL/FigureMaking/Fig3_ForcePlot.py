import numpy as np
import math
import matplotlib.pyplot as plt
import os
import thermophoreticForce as thermo
import HeightSimulation as theory

#Set directory where the file is 6-23-16-03/
os.chdir("C:/Users/Admin/Desktop/Levitating Ice Particles/code/SimulationData/")

#$$$$$ Commented out code is an alternative way to load data, also invert data to get it into poper format

#data = np.loadtxt("C:/Users/Admin/Desktop/Levitating Ice Particles/code/datatest.txt",delimiter=" ")
#hotPlateData = np.loadtxt("C:/Users/Admin/Desktop/Levitating Ice Particles/code/hotPlate.txt")
#hotPlateData = hotPlateData[::-1]
#normalizedData = data[::-1]

plateSpacing = .010
particleRadius = 13*10**-6
particleWeight = 4*math.pi/3*particleRadius**3*998*9.81
pressure = 6

x = np.linspace(0.0,1.0,600)
for file in os.listdir(os.getcwd()):
    if file.startswith("tempdata") and file.endswith(".txt"):
        normalizedData = np.loadtxt(file)
        normalizedData = normalizedData[::-1]
       # normalizedData = theory.getData(file)
        normalizedData = theory.cleanData(normalizedData)
        constants = theory.getConstants(file)
        pressure = 6
        #pressure = constants[3]
        y = thermo.calculateForce(normalizedData,particleRadius,plateSpacing,x,pressure)
       # plt.plot(x,y)
#x = np.linspace(0.0,1.0,600)
y = thermo.calculateForce(normalizedData,particleRadius,plateSpacing,x,pressure)

#weights = np.full((len(x)),particleWeight)
weights = np.full((len(x)),1.0)
#plt.plot(x*10,weights,c='k',linestyle='--')
plt.plot(x*10,thermo.temperature(normalizedData,x),c='k')
plt.show()
