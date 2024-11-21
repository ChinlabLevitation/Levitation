import numpy as np


def r_partial(values, dr):
    R_diffs = np.zeros(np.shape(values))
    R_diffs[0, :] = (values[1, :] - values[0, :])/dr
    R_diffs[-1, :] = (values[-1, :] - values[-2, :])/dr
    for j in range(1, np.shape(values)[0] - 1):
        R_diffs[j, :] = (values[j + 1, :] - values[j - 1, :]) / (2*dr)
    return R_diffs


def z_partial(values, dz):
    Z_diffs = np.zeros(np.shape(values))
    Z_diffs[:, 0] = (values[:, 1] - values[:, 0])/dz
    Z_diffs[:, -1] = (values[:, -1] - values[:, -2])/dz
    for j in range(1, np.shape(values)[1] - 1):
        Z_diffs[:, j] = (values[:, j + 1] - values[:, j - 1]) / (2*dz)
    return Z_diffs


def gradient(values, dr, dz):
    R_diffs = r_partial(values=values, dr=dr)
    Z_diffs = z_partial(values=values, dz=dz)
    return np.array([R_diffs, Z_diffs])

