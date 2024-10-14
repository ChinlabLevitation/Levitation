import numpy as np
import matplotlib.pyplot as plt
import ftSimulation as ftSim
import math
import matplotlib as mpl

mpl.rc('font',family='Times New Roman')
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20

plt.subplot(111)

path = "C:/Users/Real/Desktop/PaperData/Data/"

kn1 = np.loadtxt(path+"8-5-16_ft.txt",usecols=(0,),skiprows=1)
ft1 = np.loadtxt(path+"8-5-16_ft.txt",usecols=(1,),skiprows=1)

index = 33
tempHoldkn1 = kn1[index:]
tempHoldft1 = ft1[index:]

kn2 = np.loadtxt(path+"6-23-16-3_ft.txt",usecols=(0,),skiprows=1)
ft2 = np.loadtxt(path+"6-23-16-3_ft.txt",usecols=(1,),skiprows=1)

tempHoldkn2 = kn2[56:]
tempHoldft2 = ft2[56:]

knudsen = np.linspace(.5,1.25,600)
ftTheory = ftSim.fCompute(knudsen)

yamTheory = ftSim.yamFit(knudsen)



#freelimit =  16*math.sqrt(math.pi)/15.0
#x = np.linspace(.5,2,500)

#plt.plot(x,np.full_like(x,freelimit),linestyle='--',c='black')
plt.scatter(tempHoldkn1,tempHoldft1,label='radius = 11.3 micron',c='black',s=16,marker='o',facecolor='none')
plt.scatter(tempHoldkn2,tempHoldft2,label='radius = 13.3 micron',c='red',s=16,marker='o')
plt.plot(knudsen,ftTheory,label='New Takata Theory',c='black',linewidth=3,alpha=1)
plt.plot(knudsen,yamTheory,c='black',linewidth=3,linestyle='--')
#plt.legend(loc='upper left')
plt.xlabel(r'Knudsen number  $\mathit{Kn}$', fontsize = 21)
plt.ylabel('Hydrodynamic Factor  $\mathit{f_{T}}$', fontsize = 21)
plt.axis([.62,1.23,.2,1])

plt.yticks([.2,.4,.6,.8,1.0])
#plt.xticks([.7,.9,1.1])

plt.subplots_adjust(bottom=.12)
plt.show()