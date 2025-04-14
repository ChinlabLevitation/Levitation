import numpy as np
import math
import matplotlib.pyplot as plt

temps = [
     0.01, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15,
    -16, -17, -18, -19, -20, -21, -22, -23, -24, -25, -26, -27, -28, -29, -30, -31, -32,
    -33, -34, -35, -36, -37, -38, -39, -40, -45, -50, -55, -60, -65, -70, -75, -80
]

temps = np.array(temps)


pressures = [
    611.657, 611.15, 562.67, 517.72, 476.06, 437.47, 401.76, 368.73, 338.19, 309.98,
    283.94, 259.90, 237.74, 217.32, 198.52, 181.22, 165.30,
    150.68, 137.25, 124.92, 113.62, 103.26, 93.77, 85.10, 77.16, 69.91, 63.29, 57.25,
    51.74, 46.73, 42.16, 38.01, 34.24, 30.82,
    27.71, 24.90, 22.35, 20.04, 17.96, 16.07, 14.37, 12.84, 7.202, 3.936, 2.093, 1.080,
    0.540, 0.261, 0.122, 0.055
]

pressures = np.array(pressures)

k_B = 1.38e-23

k_B = 1.38e-23
m = 2.99e-26
T_ref = 300

pressures = pressures / math.sqrt(2 * math.pi * m * k_B * T_ref)


plt.figure(figsize=(8, 6))
plt.plot(temps, pressures, marker='o')
plt.xlabel('Temperature (°C)')
plt.ylabel('Normalized Pressure (Pa·s/m)')
plt.title('Normalized Pressure vs Temperature')
plt.grid(True)
plt.tight_layout()
plt.show()


