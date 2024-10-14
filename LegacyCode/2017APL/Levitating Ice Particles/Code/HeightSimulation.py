import numpy as np
import math
import thermophoreticForce as thermo

#$$$ This package is used to get the theoretical height values for a specific chamber configuration
#$$$ look at the plot file to see how all of this works

def getData(path):
    data = np.genfromtxt(path,skip_header=1,max_rows=1)
    data = data[::-1]
    data = cleanData(data)
    return data
    
def cleanData(data):
    data = np.array(data)
    bottomTemp = data[0]
    topTemp = data[len(data)-1]
    
    bottomTempIndices = np.where(data==bottomTemp)[0]
    topTempIndices = np.where(data==topTemp)[0]
    
    bothIndices = np.union1d(bottomTempIndices,topTempIndices)
    data = np.delete(data,bothIndices)
    data = np.insert(data,0,bottomTemp)
    data = np.insert(data,len(data),topTemp)
    return data
    
def getConstants(path):
    dataSet = np.genfromtxt(path,max_rows=1)
    return dataSet

def findTheoryHeight(data,particleRadius,pressure,plateSpacing):
    particleWeight = 4*math.pi/3*particleRadius**3*998*9.81
    x=np.linspace(0,1,200)
    force = thermo.calculateForce(data,particleRadius,plateSpacing,x,pressure)
    arrayIndex = findArrayIndex(force,particleWeight)
    theoryHeight = (arrayIndex+1)/float(200.0)
    #print(force[arrayIndex+1]-particleWeight)
    return theoryHeight
    
def findArrayIndex(arr,value):
    arr = np.array(arr)
    subtractedArray = arr-float(value)
    indices = np.where(subtractedArray>=0.0)[0]
    if len(indices)==0:
        return -1
    return np.amax(indices)
    #idx = np.searchsorted(array, value, side="left")
    #print(idx)
    #if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
    #    return idx-1
    #else:
    #    return idx  
        
def getTheoryHeightArray(path,particleRadius):
    data = getData(path)
    constants = getConstants(path)
    
    plateSpacing = constants[2]/float(1000.0)
    pressure = constants[3]
    experimentalHeight = constants[4]
    
    ###Important pressure rescaling, change when frankie updates data format
    if (experimentalHeight < 1) and (experimentalHeight > 0):
        pressure = thermo.pressureRescaling(pressure,data,experimentalHeight)
             
    topPlateTemp = data[len(data)-1]
    
    theoryHeight = findTheoryHeight(data,particleRadius,pressure,plateSpacing)   
    resultingData = np.array([topPlateTemp,experimentalHeight,theoryHeight])
    return resultingData