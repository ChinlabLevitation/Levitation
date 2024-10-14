import numpy as np
import math

def rphiz_to_xyzprime(r, phi, z, gam, eps):
    # r, phi, z are the cylindrical coordinates of the system
    # gam and eps are the rotation angles 
    # xp, yp, zp are the coordinates along axes of the cylinder
    a = np.array([[np.cos(gam)*np.cos(eps), -np.cos(gam)*np.sin(eps), np.sin(gam)], 
                  [np.sin(eps), np.cos(eps), 0], 
                  [-np.sin(gam)*np.cos(eps), np.sin(gam)*np.sin(eps), np.cos(gam)]])
    b = np.array([r, phi, z])
    [xp, yp, zp] = np.matmul(a, b)
    return xp, yp, zp # matrix multiplacation using r,phi, and z, giving us xprime,yprime,and zprime.

def xyzprime_to_rphiz(xp, yp, zp, gam,eps):
    # xp, yp, zp are the coordinates along axes of the cylinder
    # gam and eps are the rotation angles 
    # r, phi, z are the cylindrical coordinates of the system
    a = np.array([[np.cos(gam)*np.cos(eps), np.sin(eps), -np.sin(gam)*np.cos(eps)],
                  [-np.cos(gam)*np.sin(eps), np.cos(eps), np.sin(gam)*np.sin(eps)], 
                  [np.sin(gam), 0 , np.cos(gam)]])
    b = np.array([xp, yp, zp])
    [r, phi, z] = np.matmul(a, b)
    return r, phi, z # matrix multiplacation using xprime... and outputting r,phi,z

dT_dr = []
dT_dphi = 0
dT_dz = []


def Force_gradientT_xyzp(m,N,h,vis,L,R,a,T,dT_dy,dT_dz):

    a = ((-m*N*np.pi*R**2)/np.sqrt(h))*((vis)*3/4)
    b = (((8/3)+(2*L/R)*(1-(a/4))*(1/T)*(dT_dy))+(((8/3) + a*(L/R))*((1/T)*dT_dz)))
    c = a*b


