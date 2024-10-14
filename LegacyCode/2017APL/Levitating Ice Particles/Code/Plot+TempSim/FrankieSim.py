"""
1. Open .csv files to get temperature distribution data.
2. Compute the thermophoretic force at each cell along the centered Z-axis.
"""

import math
import time
import csv
#import openpyxl
import os
start_time = time.time()
init_time = time.time()

#Radius of particle in meters
radius=13*10**-6
#Weight of particle
weight=4*math.pi/3*radius**3*995*9.81
#Thermoconductivity
def k(T):
	return 0.00119+(8.94685*10**-5)*T-(2.18293*10**-8)*T**2
#Mean free path
def mfp(k,p,T):
	return 4*k/(5*p*133.322)*math.sqrt(4.81*10**-26*T/(2*1.38065*10**-23))
#Knudsen number
def Kn(mfp,radius):
	return mfp/radius
#Pressure at particle
def press_particle(p,T):
	return p*math.sqrt(T/295.0944)
#Thermophoretic force constant
def f_T(Kn):
	return (9*Kn**3)/(1+4.4844*Kn**2)/(1+Kn)
#Thermophoretic force
def force(f_T,radius,k,T,grad_T):
	return f_T*radius**2*k*grad_T/math.sqrt(2*1.38*10**-23*T/(4.8*10**-26))

#Looks into a temperature distribution and finds the theoretical height
def process(dist):
	assert isinstance(dist,list)
	#Obtain the time of the distribution
	time=dist[0][0]
	N=int(dist[0][1])
	gap=float(dist[0][2])
	press_room=float(dist[0][3])
	print "Time: "+ str(time)
	dist.pop(0)
	#Obtains the list of temperature values along the centered Z-axis.
	ls=[]
	for element in dist[0]:
		ls.append(float(element))
	result = []
	#Compute the gradient at every cell.
	l=len(ls)
	print weight
	for i in range(l):
		if i == 0:
			grad=(ls[1]-ls[0])/(gap/N)
		elif i == l-1:
			grad=(ls[l-1]-ls[l-2])/(gap/N)
		else:
			grad=(-0.5*ls[i-1]+0.5*ls[i+1])/(gap/N)**2
		f = force(f_T(Kn(mfp(k(ls[i]), press_particle(press_room, ls[i]),
							 ls[i]), radius)), radius, k(ls[i]), ls[i], grad)
		#Compiles the cell temperature, gradient and thermophoretic force
		result.append([ls[i],grad,f])
	return result
for file in os.listdir(os.getcwd()):
    if file.startswith("sim") and file.endswith(".csv"):
		print str(file)
		with open(file, 'r') as csvfile:
			dist=list(csv.reader(csvfile, delimiter=' '))
			for row in process(dist):
				print row

elapsed_time = time.time() - init_time
print
print "Total time taken: "+str(elapsed_time)+"s"