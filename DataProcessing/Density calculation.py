# Chin Lab Levitation :
# This code is for measuring the 2D surface density of 
# a material given thier dx, dy, and their density (g/m^3)


import numpy as np
import math
import matplotlib.pylab as plt
import matplotlib as mpl
from tqdm import tqdm

# number for our different materials 
# everything is in m
synth1_R = 2.0e-5/2
synth2_R = 4.3e-5/2
synth3_R = 1.9e-5/2
synth4_R = 2.7e-5/2

pi = np.pi
ro_silk = 1.3e6 # in g/m^3
ro_poly = 1.37e6 # in g/m^3

# each dx and dy length and averaging the multiple values in m

syn_Dx = [[.000172, .000172, .000169, .000170, .000176, .000168], [.000183, .000202, .000200, .000194, .000189, .000189],[.000177, .000176, .000177, .000145, .000147, .000153], [.000173, .000173, .000165, .000164, .000165, ]]
syn_Dy = [[.000138, .000133, .000136, .000112, .000114, .000114], [.000164, .000167, .000155, .000170, .000174, .000159],[.000113, .000112, .000114, .000111, .000112, .000120], [.000160, .000165, .000157, .000157, .000165]]

mean_x = [np.mean(Dx) for Dx in syn_Dx]
mean_y = [np.mean(Dy) for Dy in syn_Dy]

syn1Dx = mean_x[0]
syn1Dy = mean_y[0]

syn2Dx = mean_x[1]
syn2Dy = mean_y[1]

syn3Dx = mean_x[2]
syn3Dy = mean_y[2]

syn4Dx = mean_x[3]
syn4Dy = mean_y[3]

# fabric density equation 

if syn1Dx and syn1Dy:   #makes sure both mean_x[0] and mean_y[0] are not empty
    density1 = ro_poly * pi * (synth1_R**2) * (1/syn1Dx + 1/syn1Dy)
    print(density1)

if syn2Dx and syn2Dy: 
    density2 = ro_silk * pi * (synth2_R**2) * (1/syn2Dx + 1/syn2Dy)
    print(density2)

if syn3Dx and syn3Dy:
    density3 = ro_poly * pi * (synth3_R**2) * (1/syn3Dx + 1/syn3Dy)
    print(density3)

if syn4Dx and syn4Dy:
    density4 = ro_poly * pi * (synth4_R**2) * (1/syn4Dx + 1/syn4Dy)
    print(density4)

