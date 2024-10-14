import numpy as np
import matplotlib.pyplot as plt

path = "C:/Users/Admin/Desktop/PaperData/Data/"
temp = np.loadtxt(path+'temperaturedata.txt')
temp = temp[::,1:91:]
temp = np.fliplr(np.rot90(temp,k=0))

#levels = np.arange(77,305,50,277)
#levels = np.array([77,100,125,150,175,200,225,250,275,290,294.9999999])
#levels = np.arange(77,277,25)
levels = np.array([77,102,127,152,177,202,227,252,277,276.99])
#wanted = np.arange(0,100,1)


X2, Y2 = np.mgrid[0:21:128j, -2.3:11.7:90j]
CS = plt.contour(X2,Y2,temp,levels=levels,colors='k')


plt.xlabel('Radial Distance  $\mathit{r}$  [mm]',fontsize=16)
plt.ylabel('Height  $\mathit{h}$  [mm]',fontsize=16)
plt.clabel(CS, inline=1, fontsize=10, inlinespacing = 0,fmt='%1.0f',manual=True)

#plt.imshow(temp)
plt.show()