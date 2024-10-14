import numpy as np
import matplotlib.pyplot as plt

a = np.fromfile("C:/Users/Admin/Desktop/Levitating Ice Particles/Code/Precise.txt",dtype=float,count=-1,sep=" ")
a = a.reshape((128,128))
plt.imshow(a)
#gx,gy = np.gradient(a)
plt.contour(a,14,colors='k')
plt.axis('off')
#plt.title('Heat Distribution of Experimental Chamber')

#x = y = np.arange(0,64,1)
#X, Y = np.meshgrid(x, y,indexing='ij')
##plt.quiver(Y,X,-1*gy,-1*gx,angles='xy',scale=400)
#
#skip = (slice(None, None, 2), slice(None, None, 2))
#
#fig, ax = plt.subplots()
#ax.quiver(Y[skip], X[skip], -gy[skip], -gx[skip], a[skip],color='k')
##ax.set(aspect=1, title='Quiver Plot')
#ax.imshow(a)
plt.show()