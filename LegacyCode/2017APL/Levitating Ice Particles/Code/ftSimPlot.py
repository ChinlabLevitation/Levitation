import ftSimulation as ftSim
import HeightSimulation as heightSim
import numpy as np
import os
import matplotlib.pyplot as plt

#6-23-16-3 r = 13.3
#8-5-16 r = 11.3

#os.chdir("C:/Users/Admin/Desktop/Levitating Ice Particles/code/SimulationData/8-5-16")
os.chdir("C:/Users/Admin/Desktop/Levitating Ice Particles/code/SimulationData")
path = "C:/Users/Admin/Desktop/8-5-16.txt"

particleRadius = 11.3*10**-6
#particleWeight = 4*math.pi/3*particleRadius**3*998*9.81
#pressure = 6
kn = np.array([])
ft = np.array([])

#for file in os.listdir(os.getcwd()):
#    if file.startswith("sim") and file.endswith(".csv"):
#        temp = ftSim.simulateFt(file,particleRadius)
#        kn = np.append(kn,temp[0])
#        ft = np.append(ft,temp[1])    

#Testing
for file in os.listdir(os.getcwd()):
    os.chdir("C:/Users/Admin/Desktop/Levitating Ice Particles/code/SimulationData/"+file)
    particleRadius = np.loadtxt("C:/Users/Admin/Desktop/Levitating Ice Particles/code/SimulationData/"+file+"/radius.txt")*10**-6
    print(particleRadius)
    for file in os.listdir(os.getcwd()):
        if file.startswith("sim") and file.endswith(".csv"):
            temp = ftSim.simulateFt(file,particleRadius)
            kn = np.append(kn,temp[0])
            ft = np.append(ft,temp[1]) 
    
    kn = np.trim_zeros(kn)
    ft = np.trim_zeros(ft)

kn = np.trim_zeros(kn)
ft = np.trim_zeros(ft)

np.savetxt(path,np.transpose([kn,ft]),fmt='%.4f')

knudsen = np.linspace(0,2,600)
ftTheory = ftSim.fCompute(knudsen)

plt.scatter(kn,ft,label='Data')
plt.plot(knudsen,ftTheory,c='red',label='Theory')
plt.legend()
plt.axis([0, 2,0, 2])
plt.show()