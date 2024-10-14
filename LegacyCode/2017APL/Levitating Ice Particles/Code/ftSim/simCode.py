import math
import numpy as np
import matplotlib.pyplot as plt

def funct(kn,Aw,Ao,Hw,Ho,k21):
    x = Aw*Ho-Ao*(Hw+float(5)*math.sqrt(3.1415926)/float(4.0)*kn/k21)
    y=(Hw+float(5)*math.sqrt(3.1415926)/float(4.0)*kn/k21)**(-1)
    return 16*3.1415/float(5.0)*x*y

def f_T(knudsenNumber):
    return (9*knudsenNumber**3)/(1+4.4844*knudsenNumber**2)/float((1+knudsenNumber))

def fCompute(kn):
    result = np.empty([len(kn)])
    for n in range(len(kn)):
        result[n] = f_T(kn[n])
    return result

takata = np.loadtxt('C:/Users/Admin/Desktop/Levitating Ice Particles/Code/ftSim/takata.txt')
data = np.loadtxt('C:/Users/Admin/Desktop/Levitating Ice Particles/Code/ftSim/yamamoto.txt')


knudsen = np.empty([len(data)])
ft = np.empty([len(data)])
for n in range(0,len(data)):
    knudsen[n] = data[n,0]
    ft[n] = funct(data[n,0],data[n,1],data[n,2],data[n,3],data[n,4],.03)
    
x = np.linspace(0,50,600)
fOld = fCompute(x)

takatakn = takata[:,0]
takataft = takata[:,1]

test = 0
polyfitYam= np.poly1d(np.polyfit(knudsen[test:8],ft[test:8],2))

polyfit = np.poly1d(np.polyfit(takatakn[4:8],takataft[4:8],2))
xp = np.linspace(0,5,500)

plt.scatter(knudsen,ft,label='Yamamoto, k21=.05')
plt.scatter(takatakn,takataft,c='red',label='TakataData')
plt.plot(xp,polyfit(xp),label='Takata Fit New',c='green')
plt.plot(x,fOld,c='red',label='Takata Fit Old')
plt.plot(xp,polyfitYam(xp),c='orange',label='Yamamoto Fit')
plt.axis([0,2.5,0,1.5])
plt.legend(loc='upper left')
plt.show()

#plt.axis([0,3,0,2])