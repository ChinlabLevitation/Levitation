import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from scipy.interpolate import interp1d

mpl.rcParams.update(mpl.rcParamsDefault)
from tqdm import tqdm           # optional; a nice way to show a progress bar for long loops


# Load data
df = pd.read_csv('Ice Levitation Data(Line Data - Nathan).csv')
time_stamp = df.iloc[74:88, 2].to_numpy()
pressure_data = df.iloc[74:88, 3].astype(float).to_numpy()
bottom_plate_temperature_data = df.iloc[74:88, 4].astype(float).to_numpy() + 273.15
height_ice_data = df.iloc[74:88, 26].astype(float).to_numpy()/100
height_data = df.iloc[74:88, 15].astype(float).to_numpy()/100
height_err = 0.0005



target = 'temperature'          # 'height' or 'temperature'
target_height = 0.007           # target height from bottom plate (meters)
target_temp = 190               # target temperature (Kelvin)

iterate_through_temps = True
number_of_temps = 7
start_temp = 293
end_temp = 373

use_data = True


# temperatures in Kelvin
top_plate_temp = 77             # top bucket (liquid nitrogen)
bottom_plate_temp = 300         # bottom plate (hot plate)
boundary_temp = 290             # boundary temperature (room temperature)


show_temp_distribution = False
show_temp_slice = True




if iterate_through_temps:
    if use_data:
        bottom_plate_temps = np.linspace(np.min(bottom_plate_temperature_data), np.max(bottom_plate_temperature_data), number_of_temps)
    else:
        bottom_plate_temps = np.linspace(start_temp, end_temp, number_of_temps)
else:
    bottom_plate_temps = np.array([bottom_plate_temp])


## define physical constants and simulation parameters ##

# thermal conductivities in W/m*K
kcopper = 385
ksteel = 50.2
kairpower = 0.91823
kairconstant = 0.00014


def kair(T_):
    # this is a rough approximation of the thermal conductivity of air as a function of temperature
    # in the range we care about. See H20/Theory/Old Theory/Graph8.png. We should double-check this.
    return kairconstant * T_ ** kairpower


# geometry of the system in meters
bucket_radius = 0.017               # radius of the top plate
plate_radius = 0.02                 # radius of the bottom plate
bottom_plate_height = 0.01          # height of the bottom plate
separation_height = 0.01            # separation between the top bucket and the bottom plate
simulation_radius = 0.03            # radius of the simulation
simulation_height = 0.03            # height of the simulation


# define the number of grid points
Nr = 300                            # number of grid points in the radial direction
Nz = 300                            # number of grid points in the axial direction


##### no need to change anything below this line #####


# compute the grid spacing
dr = simulation_radius / (Nr - 1)     # grid spacing in the radial direction
dz = simulation_height / (Nz - 1)     # grid spacing in the axial direction

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
plate_mask = (R < plate_radius) & (Z < bottom_plate_height)
bucket_mask = (R < bucket_radius) & (Z > separation_height + bottom_plate_height)
air_mask = ~(plate_mask | bucket_mask)
boundary_mask = np.isclose(R, r[-1]) | (np.isclose(Z, z[0]) & (~plate_mask)) | (np.isclose(Z, z[-1]) & (~bucket_mask))


# initialize various variables
dradd = np.zeros((Nr, Nz))
drsub = np.zeros((Nr, Nz))
dzadd = np.zeros((Nr, Nz))
dzsub = np.zeros((Nr, Nz))


# define the number of iterations to perform
iterations = 20000              # number of iterations to perform


# start the simulation
T_s = []
for temp in bottom_plate_temps:
    T = np.full((Nr, Nz), boundary_temp)

    # apply boundary conditions
    T[boundary_mask] = boundary_temp
    T[plate_mask] = temp
    T[bucket_mask] = top_plate_temp

    for i in tqdm(range(iterations)):

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
        T[plate_mask] = temp
        T[bucket_mask] = top_plate_temp


    T_s.append(T.copy())


if show_temp_distribution and not iterate_through_temps:
    fullT = np.concatenate((np.flip(T_s[0][1:, :], 0), T_s[0]), axis=0)
    x = np.concatenate((-np.flip(r[1:], 0), r), axis=0)

    plt.figure(figsize=(10,4))
    cs = plt.pcolormesh(x, z, np.transpose(fullT), cmap='plasma', vmin=top_plate_temp, vmax=bottom_plate_temps[0])

    plt.colorbar(label='Temperature $T$ (K)', ticks=[100, 150, 200, 250, 300, 350])
    plt.gca().set_aspect('equal')

    cs2 = plt.contour(x, z, np.transpose(fullT), levels=np.linspace(top_plate_temp, bottom_plate_temps[0], 15), colors='black',linewidths=1)

    plt.xlabel('Position $x$ (m)')
    plt.ylabel('Height $z$ (m)')
    plt.show()

# Extract the temperature slice at r=0 for the selected z-range
z_indices = np.where((z >= bottom_plate_height) & (z <= separation_height + bottom_plate_height))[0]
center_slices = [T_[0, z_indices] for T_ in T_s]
z_slice = z[z_indices] - bottom_plate_height



if show_temp_slice and not iterate_through_temps and target == 'height':
    plt.figure(figsize=(8, 5))
    plt.plot(center_slices[0], z_slice, 'b-', label='Temperature at r=0')
    plt.xlabel('Temperature $T$ (K)')
    plt.ylabel('Height $z$ (m)')
    plt.title('Height vs Temperature at Center (r=0) Between Plates')
    plt.grid(True)
    plt.legend()
    plt.show()

# Plot the 1D temperature slice
if show_temp_slice and not iterate_through_temps and target == 'temperature':
    plt.figure(figsize=(8, 5))
    plt.plot(z_slice, center_slices[0], 'b-', label='Temperature at r=0')
    plt.xlabel('Height $z$ (m)')
    plt.ylabel('Temperature $T$ (K)')
    plt.title('Temperature vs Height at Center (r=0) Between Plates')
    plt.grid(True)
    plt.legend()
    plt.show()



heights_at_temp = []
temps_at_height = []

for temp_slice in center_slices:
    if target == 'temperature':
        order = np.argsort(temp_slice)
        interp_func = interp1d(temp_slice[order], z_slice[order], kind='cubic',bounds_error=False, fill_value="extrapolate")
        heights_at_temp.append(interp_func(target_temp))
    else:
        order = np.argsort(z_slice)
        interp_func = interp1d(z_slice[order], temp_slice[order], kind='cubic', bounds_error=False,fill_value="extrapolate")
        temps_at_height.append(interp_func(target_height))


#Plot temperature vs bottom plate temperature of an isotherm
if show_temp_slice and iterate_through_temps and target == 'temperature':
    plt.figure(figsize=(8, 5))
    plt.plot(bottom_plate_temps, np.array(heights_at_temp), 'b-', label=f'{target_temp:.2f}K Isotherm')
    if use_data:
        plt.errorbar(bottom_plate_temperature_data, height_data, xerr=0, yerr=height_err, fmt='o', capsize=3, label='Data Points')
    plt.xlabel('Bottom Plate Temperature (K)')
    plt.ylabel('Height $z$ (m)')
    plt.title('Height vs Bottom Plate Temperature of Isotherm')
    plt.grid(True)
    plt.legend()
    plt.show()


#Plot height vs bottom plate temperature of an isoheight
if show_temp_slice and iterate_through_temps and target == 'height':
    plt.figure(figsize=(8, 5))
    plt.plot(bottom_plate_temps, np.array(temps_at_height), 'b-', label=f'{target_height:.2f}m Isoheight')
    plt.xlabel('Bottom Plate Temperature (K)')
    plt.ylabel('Temperature $T$ (K)')
    plt.title('Temperature vs Bottom Plate Temperature of Isotherm')
    plt.grid(True)
    plt.legend()
    plt.show()

