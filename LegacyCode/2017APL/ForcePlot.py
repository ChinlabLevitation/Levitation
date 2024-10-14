import numpy as np
import math
import matplotlib.pyplot as plt
import os
import thermophoreticForce as thermo
import HeightSimulation as theory

#Set directory where the file is
os.chdir("C:/Users/Admin/Desktop/Levitating Ice Particles/code/SimulationData/6-23-16-03/")

#$$$$$ Commented out code is an alternative way to load data, also invert data to get it into poper format

#data = np.loadtxt("C:/Users/Admin/Desktop/Levitating Ice Particles/code/datatest.txt",delimiter=" ")
#hotPlateData = np.loadtxt("C:/Users/Admin/Desktop/Levitating Ice Particles/code/hotPlate.txt")
#hotPlateData = hotPlateData[::-1]
#normalizedData = data[::-1]

plateSpacing = .01036
particleRadius = 13.3*10**-6
particleWeight = 4*math.pi/3*particleRadius**3*998*9.81
pressure = 3.5

x = np.linspace(0.0,1.0,600)
for file in os.listdir(os.getcwd()):
    if file.startswith("sim") and file.endswith(".csv"):
        normalizedData = theory.getData(file)
        normalizedData = theory.cleanData(normalizedData)
        constants = theory.getConstants(file)
        pressure = constants[3]
        y = thermo.calculateForce(normalizedData,particleRadius,plateSpacing,x,pressure)
        plt.plot(x,y)
#x = np.linspace(0.0,1.0,600)
y = thermo.calculateForce(normalizedData,particleRadius,plateSpacing,x,pressure)

weights = np.full((len(x)),particleWeight)
plt.plot(x,weights)
plt.plot(x,y)
plt.show()
