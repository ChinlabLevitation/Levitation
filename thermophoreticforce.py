import numpy as np

# This file contains functions and constants used to calculate the thermophoretic force on a particle in a gas.
# Mostly, the contents of this file are copied from the 2017APL code, with some minor modifications.
# To be honest, mostly I'm just adding documentation.


# average molecular mass of air in kg, from Wolfram Alpha
molecMassAir = 4.809e-26 

# Boltzmann constant in J/K
boltzmannConstant = 1.38065e-23 

# Thermal conductivity of air in W/(m*K)
def airThermalConductivity(temp):
    return 0.00119+(8.94685*10**-5)*temp-(2.18293*10**-8)*(temp)**2

def pressure(press,temp):
	#return press*math.sqrt(temp/295.0944)
	return press

def pressureRescaling(press,data,height):
    return press*math.sqrt(295.0944/temperature(data,height))

def meanFreePath(airThermalConductivity,pressure,temp):
    return 4*airThermalConductivity/(5*pressure*133.322)*math.sqrt(molecMassAir*temp/(2*boltzmannConstant))
    
def knudsenNumber(meanFreePath,particleRadius):
    return meanFreePath/float(particleRadius)

def computeKnudsenNumber(data,height,expPress,radius):
    if (height>1) or (height < 0):
        return 0.0
    temp = temperature(data,height)
    aTC = airThermalConductivity(temp)
    press = pressure(expPress,temp)
    print(press)
    mFP = meanFreePath(aTC,press,temp)
    return knudsenNumber(mFP,radius)

def f_T(knudsenNumber):
    return (9*knudsenNumber**3)/(1+4.4844*knudsenNumber**2)/float((1+knudsenNumber))

def newf_T(knudsenNumber):
    ft = -.2906*knudsenNumber**2 + 1.354*knudsenNumber - .2706
    if (knudsenNumber < .35) or (knudsenNumber > 2):
        print("Warning: f_T inaccurate, switch from newf_T function.")
        return ft
    return ft

def thermophoreticForce(f_T,particleRadius,airThermalConductivity,temp,tempGradient):
    return -(f_T*particleRadius**2*airThermalConductivity* \
        tempGradient/math.sqrt(2*boltzmannConstant*temp/float(molecMassAir)))

#Calculates thermophoretic force for numpy array of heights, returns numpy array of forces    
def calculateForce(data,particleRadius,plateSpacing,heights,press):
    computations = len(heights)
    forceValues = np.zeros(computations)
    
    for n in range(0,computations):
        temp = temperature(data,heights[n])
        p = pressure(press,temp)
        aTC = airThermalConductivity(temp)
        mFP = meanFreePath(aTC,p,temp)
        kn = knudsenNumber(mFP,particleRadius)
        ft = f_T(kn)
        tempGrad = tempGradient(data,heights[n],plateSpacing)
        
        forceValues[n]=thermophoreticForce(ft,particleRadius,aTC,temp,tempGrad)
    return forceValues

#Calculates gradient force for numpy array of heights, returns numpy array of gradients
def calculateGradient(data,heights,plateSpacing):
    y = np.zeros(len(heights))
    for n in range(0,len(heights)):
        y[n] = tempGradient(data,heights[n],plateSpacing)
    return y   