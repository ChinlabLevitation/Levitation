import thermophoreticForce as thermo
import HeightSimulation as heightSim
import numpy as np
import math

molecMassAir = 4.81*math.pow(10,-26)
boltzmannConstant = 1.38065*math.pow(10,-23)

def calculateFt(data,particleRadius,plateSpacing,heights,press):
        heights = 1.0-heights
        temp = thermo.temperature(data,heights)
        p = thermo.pressure(press,temp)
        aTC = thermo.airThermalConductivity(temp)
        mFP = thermo.meanFreePath(aTC,p,temp)
        kn = thermo.knudsenNumber(mFP,particleRadius)
        tempGrad = thermo.tempGradient(data,heights,plateSpacing)      
        force = -(float(particleRadius)**2*aTC*tempGrad/math.sqrt(2*boltzmannConstant*temp/float(molecMassAir)))
        
        particleWeight = 4*math.pi/3*particleRadius**3*998*9.81
        
        fT = particleWeight/float(force)
        results = np.array([kn,fT])
        
        return results
def simulateFt(path,particleRadius):
    data = heightSim.getData(path)
    data = heightSim.cleanData(data)
    constants = heightSim.getConstants(path)
    
    plateSpacing = constants[2]/float(1000.0)
    pressure = constants[3]
    experimentalHeight = constants[4]
    
    if (experimentalHeight < 0) or (experimentalHeight > 1):
        return np.array([0,0])
        
    return calculateFt(data,particleRadius,plateSpacing,experimentalHeight,pressure)
    
def f_T(knudsenNumber):
    return (9*knudsenNumber**3)/(1+4.4844*knudsenNumber**2)/float((1+knudsenNumber))

def newf_T(knudsenNumber):
    ft = -.2906*knudsenNumber**2 + 1.354*knudsenNumber - .2706
    if (knudsenNumber < .35) or (knudsenNumber > 2):
        print("Warning: f_T inaccurate, switch from newf_T function.")
        return ft
    return ft

def fCompute(kn):
    result = np.empty([len(kn)])
    for n in range(len(kn)):
        result[n] = newf_T(kn[n])
    return result