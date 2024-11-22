import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
import math

molecMassAir = 4.81 * 10**-26  # kg
boltzmannConstant = 1.38065 * 10**-23  # J/K
particleRadius = 18.3 * 10**-6  # m
particleWeight = (4 * math.pi / 3) * particleRadius**3 * 998 * 9.81  # N (density * volume * g)
a = 0.00014  # W/(m·K)
b = 0.91823  # dimensionless

def airThermalConductivity(temp):
    return a * temp ** b  # W/(m·K)

def meanFreePath(pressure, temp):
    # Returns mean free path in meters
    return 4*airThermalConductivity(temp)/(5*pressure)*np.sqrt(molecMassAir*temp/(2*boltzmannConstant))

def knudsenNumber(pressure, temp, particleRadius):
    # Dimensionless number
    return meanFreePath(pressure, temp)/float(particleRadius)

def f_T(pressure, temp, particleRadius):
    # Dimensionless thermal force coefficient
    K_n = knudsenNumber(pressure, temp, particleRadius)
    return (9 * K_n ** 3) / (1 + 4.4844 * K_n ** 2) / (1 + K_n.astype(float))


def odes(z, y):
    T = y[0]
    dTdz = y[1]
    d2Tdz2 = - (b * dTdz**2) / T
    return np.vstack((dTdz, d2Tdz2))


def bc(ya, yb):
    return np.array([ya[0] - 293, yb[0] - 77])  # K


z = np.linspace(0, 0.01, 1000)  # m
dz = z[1] - z[0]  # m

T_guess = np.linspace(293, 77, z.size)  # K
dTdz_guess = np.gradient(T_guess, z)  # K/m
y_guess = np.vstack((T_guess, dTdz_guess))

solution = solve_bvp(odes, bc, z, y_guess)
T = np.array(solution.sol(z)[0])  # K
grad_T = (T[1:] - T[:-1])/dz  # K/m

# Find index closest to z = 0.5cm = 0.005m
z_index = np.abs(z - 0.005).argmin()

# Create pressure range to plot against
pressures = np.linspace(0.5, 10, 100) * 133.322  # Convert from Torr to Pascal (Pa)

# Calculate force at z=0.5cm for different pressures
forces = []
for p in pressures:
    # Force calculation in Newtons (N)
    # f_T: dimensionless
    # airThermalConductivity: W/(m·K)
    # particleRadius^2: m^2
    # sqrt(2kT/m): m/s
    # grad_T: K/m
    F = -f_T(p, T[z_index], particleRadius) * (a * T[z_index]**b) * particleRadius**2/np.sqrt(2 * boltzmannConstant * T[z_index]/molecMassAir) * grad_T[z_index]
    forces.append(F)

forces = np.array(forces)  # N
mg = particleWeight  # N


from tqdm import tqdm  # optional; a nice way to show a progress bar for long loops

# thermal conductivities in W/m*K
kcopper = 385
ksteel = 50.2
kairpower = 0.918
kairconstant = 0.00014


def kair(T):
    # this is a rough approximation of the thermal conductivity of air as a function of temperature
    # in the range we care about. See H20/Theory/Old Theory/Graph8.png. We should double check this.
    return kairconstant * T ** kairpower


# temperatures in K
bucket_temp = 73  # top bucket (liquid nitrogen)
plate_temp = 315.5  # bottom plate (hot plate)
boundary_temp = 289  # boundary temperature (room temperature)

# geometry of the system in inches
bucket_radius = .367  # radius of the top bucket
plate_radius = .665  # radius of the bottom plate
plate_height = 1 - 0.6102362  # height of the bottom plate
separation_height = 0.6102362 - 0.1969  # separation between the top bucket and the bottom plate
simulation_radius = 1  # radius of the simulation
simulation_height = 1  # height of the simulation

# define the number of grid points
Nr = 128  # number of grid points in the radial direction
Nz = 128  # number of grid points in the axial direction

## no need to change anything below this line ##

# convert geometry to meters
bucket_radius *= 0.0254
plate_radius *= 0.0254
plate_height *= 0.0254
separation_height *= 0.0254
simulation_radius *= 0.0254
simulation_height *= 0.0254

# compute the grid spacing
dr = simulation_radius / Nr  # grid spacing in the radial direction
dz = simulation_height / Nz  # grid spacing in the axial direction

# Create 2D grid for z and r coordinates
z = np.linspace(0, simulation_height, Nz)  # Height in meters
r = np.linspace(0, simulation_radius, Nr)  # Radius in meters
Z, R = np.meshgrid(z, r)

# Initialize arrays
T = np.full_like(Z, boundary_temp)  # Initialize all points to boundary temp
k = np.zeros_like(Z)

# Create masks for different regions
plate_mask = np.logical_and(R < plate_radius, Z < plate_height)
bucket_mask = np.logical_and(R < bucket_radius, Z > separation_height + plate_height)
air_mask = np.logical_not(np.logical_or(plate_mask, bucket_mask))
boundary_mask = R >= simulation_radius - dr  # Only side walls are boundary conditions

# set up the initial thermal conductivity
k[plate_mask] = kcopper
k[bucket_mask] = ksteel
k[air_mask] = kair(T[air_mask])

# initialize various variables
dradd = np.zeros((Nr, Nz))
drsub = np.zeros((Nr, Nz))
dzadd = np.zeros((Nr, Nz))
dzsub = np.zeros((Nr, Nz))

# define the number of iterations to perform
iterations = 10000  # number of iterations to perform

# start the simulation
for t in tqdm(range(iterations)):
    # apply boundary conditions
    T[boundary_mask] = boundary_temp
    T[plate_mask] = plate_temp
    T[bucket_mask] = bucket_temp

    # compute thermal conductivity
    k = kair(T)

    # compute terms in the bulk
    dradd = k / (dr ** 2)
    dzadd = k / (dz ** 2)
    drsub[1:-1, :] = k[1:-1, :] / (2 * dr * R[1:-1, :]) + (k[2:, :] - k[:-2, :]) / (2 * dr) ** 2
    dzsub[:, 1:-1] = (k[:, 2:] - k[:, :-2]) / (2 * dz) ** 2

    # compute terms at the boundaries
    drsub[0, :] = k[0, :] / (dr ** 2)
    drsub[-1, :] = - k[-1, :] / (dr ** 2)
    dzsub[:, 0] = k[:, 0] / (dz ** 2)
    dzsub[:, -1] = - k[:, -1] / (dz ** 2)

    # compute update
    T = (np.roll(T, -1, 0) * (dradd + drsub) + np.roll(T, 1, 0) * (dradd - drsub)
         + np.roll(T, -1, 1) * (dzadd + dzsub) + np.roll(T, 1, 1) * (dzadd - dzsub))
    T = T / (2 * k * (1 / dr ** 2 + 1 / dz ** 2))

# Create symmetric temperature distribution for plotting
fullT = np.concatenate((np.flip(T[1:, :], 0), T), axis=0)
x = np.concatenate((-np.flip(r[1:], 0), r), axis=0)

plt.figure(figsize=(10, 4))
cs = plt.pcolormesh(x, z, np.transpose(fullT), cmap='plasma', vmin=bucket_temp, vmax=plate_temp)

plt.colorbar(label='Temperature $T$ (K)', ticks=[100, 150, 200, 250, 300, 350])
plt.gca().set_aspect('equal')

cs2 = plt.contour(x, z, np.transpose(fullT), levels=np.linspace(bucket_temp, plate_temp - .01, 16), colors='black', linewidths=1)

plt.xlabel('Position $x$ (m)')
plt.ylabel('Height $z$ (m)')
plt.show()



# Calculate temperature gradients
dTdr = np.zeros_like(T)
dTdz = np.zeros_like(T)

dTdr[1:-1,:] = (T[2:,:] - T[:-2,:])/(2*dr)
dTdz[:,1:-1] = (T[:,2:] - T[:,:-2])/(2*dz)

# Calculate temperature gradients for full domain
fullDTdr = np.concatenate((-np.flip(dTdr[1:,:], 0), dTdr), axis=0)
fullDTdz = np.concatenate((np.flip(dTdz[1:,:], 0), dTdz), axis=0)
gradT_mag = np.sqrt(fullDTdr**2 + fullDTdz**2)

plt.figure(figsize=(10, 4))
plt.pcolormesh(x, z, np.transpose(fullDTdz), cmap='viridis')
plt.colorbar(label='fullDTdz')
plt.xlabel('Position $x$ (m)')
plt.ylabel('Height $z$ (m)')
plt.title("fullDTdz")
plt.gca().set_aspect('equal')
plt.show()

'''

plt.figure(figsize=(10,4))
cs = plt.pcolormesh(x, z, np.transpose(gradT_mag), cmap='plasma')
plt.colorbar(label='Temperature Gradient Magnitude (K/m)')
plt.gca().set_aspect('equal')

# Add contour lines
cs2 = plt.contour(x, z, np.transpose(gradT_mag), levels=np.linspace(0, np.max(gradT_mag), 10), colors='black', linewidths=1)

plt.xlabel('Position $x$ (m)')
plt.ylabel('Height $z$ (m)')
plt.title('Temperature Gradient Field')
plt.show()



# Calculate temperature gradients
dTdr = np.zeros_like(T)
dTdz = np.zeros_like(T)

dTdr[1:-1,:] = (T[2:,:] - T[:-2,:])/(2*dr)
dTdz[:,1:-1] = (T[:,2:] - T[:,:-2])/(2*dz)

# Calculate boundary gradients
dTdr[0,:] = (T[1,:] - T[0,:])/dr
dTdr[-1,:] = (T[-1,:] - T[-2,:])/dr
dTdz[:,0] = (T[:,1] - T[:,0])/dz
dTdz[:,-1] = (T[:,-1] - T[:,-2])/dz


# Calculate Knudsen number and thermophoretic function
averagePressure = 5 * 133.322  # Convert 5 Torr to Pascal
Kn = meanFreePath(averagePressure, T)
f_t = (9 * Kn**3)/(1 + 4.4844 * Kn**2)/(1 + Kn)

# Calculate force components
Fr = -f_t * (kair(T)) * particleRadius**2/np.sqrt(2 * boltzmannConstant * T/molecMassAir) * dTdr
Fz = -f_t * (kair(T)) * particleRadius**2/np.sqrt(2 * boltzmannConstant * T/molecMassAir) * dTdz

# Create full arrays for plotting
fullFr = np.concatenate((-np.flip(Fr[1:,:], 0), Fr), axis=0)
fullFz = np.concatenate((np.flip(Fz[1:,:], 0), Fz), axis=0)
F_mag = np.sqrt(fullFr**2 + fullFz**2)

# Plot force magnitude
plt.figure(figsize=(10,4))
cs = plt.pcolormesh(x, z, np.transpose(F_mag), cmap='viridis')
plt.colorbar(label='Force magnitude (N)')
plt.gca().set_aspect('equal')



plt.xlabel('Position $x$ (m)')
plt.ylabel('Height $z$ (m)')
plt.title('Thermophoretic Force Field')
plt.show()




# Create a new figure for force magnitude and direction plot
plt.figure(figsize=(10,4))

# Plot force magnitude as background color
cs = plt.pcolormesh(x, z, np.transpose(F_mag), cmap='viridis')
plt.colorbar(label='Force magnitude (N)')

# Plot force vectors with arrows, downsampling for clarity
skip = 5  # Plot every 5th point to avoid overcrowding
plt.quiver(x[::skip], z[::skip],
          np.transpose(fullFr)[::skip,::skip]/F_mag.T[::skip,::skip],  # Normalize to get direction
          np.transpose(fullFz)[::skip,::skip]/F_mag.T[::skip,::skip],
          F_mag.T[::skip, ::skip],  # Use F_mag to color the arrows
           scale=28,  # Adjust this value to get the desired arrow length
           width=0.005,  # Adjust width if needed
           headwidth=3,
           headlength=5,
           headaxislength=4.5,
           cmap='viridis')


plt.xlabel('Position $x$ (m)')
plt.ylabel('Height $z$ (m)')
plt.title('Thermophoretic Force Field Direction')
plt.gca().set_aspect('equal')
plt.show()



# Create a new figure for force magnitude and direction plot with varying arrow sizes
plt.figure(figsize=(10,4))

# Plot force magnitude as background color
cs = plt.pcolormesh(x, z, np.transpose(F_mag), cmap='viridis')
plt.colorbar(label='Force magnitude (N)')

plt.quiver(x[::skip], z[::skip],
          np.transpose(fullFr)[::skip,::skip]/F_mag.T[::skip,::skip],  # Normalize to get direction
          np.transpose(fullFz)[::skip,::skip]/F_mag.T[::skip,::skip],
          color='white',
           scale=28,  # Adjust this value to get the desired arrow length
           width=0.005,  # Adjust width if needed
           headwidth=3,
           headlength=5,
           headaxislength=4.5,
           cmap='viridis')

plt.xlabel('Position $x$ (m)')
plt.ylabel('Height $z$ (m)')
plt.title('Thermophoretic Force Field (Arrow size proportional to force)')
plt.gca().set_aspect('equal')
plt.show()




# bucket_radius = .367 # radius of the top bucket
# plate_radius = .665 # radius of the bottom plate
# plate_height = 1 - 0.6102362 # height of the bottom plate
# separation_height = 0.6102362 - 0.1969 # separation between the top bucket and the bottom plate
# simulation_radius = 1 # radius of the simulation
# simulation_height = 1 # height of the simulation


# Create a new figure for force magnitude and direction plot with varying arrow sizes
plt.figure(figsize=(10,4))

# Plot force magnitude as background color
cs = plt.pcolormesh(x, z, np.transpose(F_mag), cmap='viridis')
plt.colorbar(label='Force magnitude (N)')

# Draw boundary lines
plt.plot([-0.0665/4, 0.0665/4], [0.01, 0.01], color='red', linestyle='-', label='Bottom Plate Boundary')
plt.plot([-0.0367/4, 0.0367/4], [0.02, 0.02], color='red', linestyle='-', label='Top Plate Boundary')

plt.quiver(x[::skip], z[::skip],
          np.transpose(fullFr)[::skip,::skip],  # NO NORMALIZATION
          np.transpose(fullFz)[::skip,::skip],
          color='white',
           scale=8*10**-22,  # Adjust this value to get the desired arrow length
           width=0.005,  # Adjust width if needed
           headwidth=3,
           headlength=5,
           headaxislength=4.5,
           cmap='viridis')

plt.xlabel('Position $x$ (m)')
plt.ylabel('Height $z$ (m)')
plt.title('Thermophoretic Force Field (Arrow size proportional to force)')
plt.gca().set_aspect('equal')
plt.legend()
plt.show()
'''