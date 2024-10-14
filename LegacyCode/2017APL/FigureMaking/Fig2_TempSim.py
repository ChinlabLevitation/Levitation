import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pylab as pl

mpl.rc('font',family='Times New Roman')
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20

plt.subplot(111)
# Load Data
path = "C:/Users/Real/Desktop/PaperData/Data/"
temp = np.loadtxt(path+'temperaturedata.txt')

#Trim Data
temp = temp[::,1:91:]
#Create a boundary at 298K at the edge to appear in graphic
zeros = np.zeros((5,90))
zeros.fill(298.00)
#Add the boundary to the data set, rotate the data set so it appears properly
temp = np.concatenate((temp,zeros))
temp = np.fliplr(np.rot90(temp,k=0))
#Specify contours that we will draw
levels = np.array([77,102,127,152,177,202,227,252,276.99,277,297.99])
#Set the size of the x,y ticks

#Plot Specific Ticks

#Create Grid, second number provides dimension on graph
X2, Y2 = np.mgrid[0:22:133j, -2.3:11.7:90j]
#Plot the contours
CS = plt.contour(X2,Y2,temp,levels=levels,colors='k')
#Add text
plt.text(6, -1.5, '277 K', fontsize=20,color='white')
plt.text(6, 10.5, '77 K', fontsize=20,color='white')
#plt.text(17, 4, '277', fontsize=18,color='black')
#Make x values for shading and shade in the plates
x = np.arange(0,16.85,.05)
x2 = np.arange(0,12.67,.05)

plt.fill_between(x,-2.3,.05,color='grey')
plt.fill_between(x2,9.98,11.7,color='grey')

plt.subplots_adjust(bottom=.12)

#Plot labels and manually edit the curves
plt.xlabel('Radial Distance  $\mathit{r}$  [mm]',fontsize=21)
plt.ylabel('Height  $\mathit{z}$  [mm]',fontsize=21)

plt.xticks([0,10,20])
plt.yticks([0,5,10])

plt.arrow(.85,.75,0,3,fc="k",ec="k",head_width=.65,head_length=.65,lw=4)
plt.arrow(8,5.5,-4,0,fc="k",ec="k",head_width=.65,head_length=.65,lw=4)

plt.clabel(CS, inline=1, fontsize=20, inlinespacing = 0,fmt='%1.0f')

plt.show()