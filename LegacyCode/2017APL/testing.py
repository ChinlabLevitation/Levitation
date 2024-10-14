import matplotlib.pyplot as plt
import thermophoreticForce as force
import numpy as np

t = np.array([77,145,185,215,240,260,270])

x = np.linspace(0,1,600)
y=np.empty(len(x))

y = force.calculateGradient(t,x,.01)

plt.plot(x,y)
plt.show()