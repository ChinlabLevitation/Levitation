import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl
from scipy.interpolate import interp1d

mpl.rcParams.update(mpl.rcParamsDefault)
from tqdm import tqdm           # optional; a nice way to show a progress bar for long loops



height = 0.007                  # height from bottom plate (meters)

# temperatures in K
bucket_temp = 77                # top bucket (liquid nitrogen)
plate_temp = 300                # bottom plate (hot plate)
boundary_temp = 290             # boundary temperature (room temperature)

show_temp_distribution = True
show_temp_slice = True


## define physical constants and simulation parameters ##

# thermal conductivities in W/m*K
kcopper = 385
ksteel = 50.2
kairpower = 0.91823
kairconstant = 0.00014


def kair(T_):
    # this is a rough approximation of the thermal conductivity of air as a function of temperature
    # in the range we care about. See H20/Theory/Old Theory/Graph8.png. We should double check this.
    return kairconstant * T_ ** kairpower


# geometry of the system in meters
bucket_radius = 0.017           # radius of the top plate
plate_radius = 0.02             # radius of the bottom plate
bottom_plate_height = 0.01      # height of the bottom plate
separation_height = 0.01        # separation between the top bucket and the bottom plate
simulation_radius = 0.03        # radius of the simulation
simulation_height = 0.03        # height of the simulation


# define the number of grid points
Nr = 300                        # number of grid points in the radial direction
Nz = 300                        # number of grid points in the axial direction


##### no need to change anything below this line #####


# compute the grid spacing
dr = simulation_radius / Nr     # grid spacing in the radial direction
dz = simulation_height / Nz     # grid spacing in the axial direction

# set up the grid
r = np.linspace(0, simulation_radius, Nr)
z = np.linspace(0, simulation_height, Nz)
Z, R = np.meshgrid(z, r)
Zpix, Rpix = np.meshgrid(np.arange(Nz), np.arange(Nr))
# meshgrid returns two 2D arrays, R and Z, that contain the r and z values at each point in the grid


# set up the temperature and thermal conductivity arrays
T = np.zeros((Nr, Nz))          # temperature
k = np.zeros((Nr, Nz))          # thermal conductivity

# set up "masks", which are arrays that are 1 where a condition is true and 0 where it is false
# these masks will be used to apply boundary conditions and to define the regions of the system
# replaces the "thermal conductivity" and "boundary_condition" functions from the old code
plate_mask = np.logical_and(R < plate_radius, Z < bottom_plate_height)
bucket_mask = np.logical_and(R < bucket_radius, Z > separation_height + bottom_plate_height)
air_mask = np.logical_not(np.logical_or(plate_mask, bucket_mask))
boundary_mask = R == simulation_radius

# set up the initial temperature distribution
T[:, :] = boundary_temp

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
iterations = 30000              # number of iterations to perform

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

# # apply boundary conditions
T[boundary_mask] = boundary_temp
T[plate_mask] = plate_temp
T[bucket_mask] = bucket_temp


if show_temp_distribution:
    fullT = np.concatenate((np.flip(T[1:, :], 0), T), axis=0)
    x = np.concatenate((-np.flip(r[1:], 0), r), axis=0)

    plt.figure(figsize=(10,4))
    cs = plt.pcolormesh(x, z, np.transpose(fullT), cmap='plasma', vmin=bucket_temp, vmax=plate_temp)

    plt.colorbar(label='Temperature $T$ (K)', ticks=[100, 150, 200, 250, 300, 350])
    plt.gca().set_aspect('equal')

    cs2 = plt.contour(x, z, np.transpose(fullT), levels=np.linspace(bucket_temp, plate_temp, 15), colors='black',linewidths=1)

    plt.xlabel('Position $x$ (m)')
    plt.ylabel('Height $z$ (m)')
    plt.show()

# Extract the temperature slice at r=0 for the selected z-range
z_indices = np.where((z >= bottom_plate_height) & (z <= separation_height + bottom_plate_height))[0]
center_slice = T[0, z_indices]
z_slice = z[z_indices] - bottom_plate_height
interp_func = interp1d(z_slice, center_slice, kind='cubic', fill_value="extrapolate")
temp_at_height = interp_func(height)
print(temp_at_height)


# Plot the 1D temperature slice
if show_temp_slice:
    plt.figure(figsize=(8, 5))
    plt.plot(z_slice, center_slice, 'b-', label='Temperature at r=0')
    plt.xlabel('Height $z$ (m)')
    plt.ylabel('Temperature $T$ (K)')
    plt.title('Temperature Profile at Center (r=0) Between Plates')
    plt.grid(True)
    plt.legend()
    plt.show()












