import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
import math

molecMassAir = 4.81 * 10**-26
boltzmannConstant = 1.38065 * 10**-23

averagePressure = 10 * 133.322
particleRadius = 13e-6
particleDensity = 917

particleWeight = (4 * math.pi / 3) * particleRadius**3 * particleDensity * 9.81

top_temp = 77
bottom_temp = 300.65

a = 0.00014
b = 0.91823

def airThermalConductivity(temp):
    return a * temp ** b
    # return 0.00119+(8.94685*10**-5)*temp-(2.18293*10**-8)*temp**2

def meanFreePath(pressure, temp):
    return 4*airThermalConductivity(temp)/(5*pressure)*np.sqrt(molecMassAir*temp/(2*boltzmannConstant))

def knudsenNumber(pressure, temp, particleRadius):
    return meanFreePath(pressure, temp)/float(particleRadius)

def f_T(pressure, temp, particleRadius):
    K_n = knudsenNumber(pressure, temp, particleRadius)
    return (9 * K_n ** 3) / (1 + 4.4844 * K_n ** 2) / (1 + K_n.astype(float))


def odes(z, y):
    T = y[0]
    dTdz = y[1]
    d2Tdz2 = - (b * dTdz**2) / T
    return np.vstack((dTdz, d2Tdz2))


def bc(ya, yb):
    return np.array([ya[0] - bottom_temp, yb[0] - top_temp])


z = np.linspace(0, 0.01, 1000)
dz = z[1] - z[0]

T_guess = np.linspace(bottom_temp, top_temp, z.size)
dTdz_guess = np.gradient(T_guess, z)
y_guess = np.vstack((T_guess, dTdz_guess))

solution = solve_bvp(odes, bc, z, y_guess)


T = np.array(solution.sol(z)[0])

grad_T = (T[1:] - T[:-1])/dz

# plt.scatter(z[:-1], f_T(averagePressure, T[:-1], particleRadius))
# plt.show()

F = -f_T(averagePressure, T[:-1], particleRadius) * (a * T[:-1]**b) * particleRadius**2/np.sqrt(2 * boltzmannConstant * T[:-1]/molecMassAir) * grad_T
mg = np.ones(np.shape(z[:-1])) * particleWeight

# Plot multiple lines
plt.figure(figsize=(10, 6))

plt.plot(z[:-1], F, label='F')
plt.plot(z[:-1], mg, label='mg')

plt.xlabel('z')
plt.ylabel('Function values')
plt.title('Multiple Lines Example')
plt.legend()
plt.grid(True)
plt.show()




