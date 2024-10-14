"""
1. Converts temperature and pressure data into an excel spreadsheet.
2. Converts tracking data into an excel spreadsheet.
3. Runs simulation and outputs results in a new folder.
4. Outputs a table of T vs grad T
"""

import math
import time
import csv
import openpyxl
import os
start_time = time.time()
init_time = time.time()

#Parameters for setup
fac = 0.6/12.7	#Scaling factor
r1 = 25.4/2		#Top plate radius
r2 = 33.7/2		#Bottom plate radius
gap = 10.00		#Gap size
N = 128			#Grid resolution of simulation
M = 100			#Number of simulations
T1 = 77			#Temp of top plate
T2 = 300		#Temp of bottom plate
T3 = 300		#Temp of surroundings
errTop = 0.01	#Maximum error; the smaller the slower the simulation
r1 = r1*fac		#N2 plate radius
z1 = 0.094		#Top plate
r2 = r2*fac		#Cu plate radius
z2 = z1+gap*fac	#Bottom plate height
err = 100
alpha = 0.918

#Paramaters for video
fr = 15 		#Frame rate
vidz = 1280		#Video size in z direction
vidx = 720		#Video size in x direction
topz = 520		#Position of top plate
botz = 520+326		#Position of bottom plate

#Functions for converting from voltage data
def pressure(x):
	return (-4.62+4.6742*10**-34*math.exp(26.1567*x))
def thermistor(x):
	return	(0.003357042+0.00025214848*math.log(x)+0.0000033743283*math.log(x)**2
			-0.000000064957311*math.log(x)**3)**-1
def thermocouple(x):
	return	(2.5173462*10**-2*(x/480*10**6)-1.1662878*10**-6*(x/480*10**6)**2
			-1.083363*10**-9*(x/480*10**6)**3-8.977354*10**-13*(x/480*10**6)**4
			-3.7342377*10**-16*(x/480*10**6)**5-8.6632643*10**-20*(x/480*10**6)**6
			-1.0450598*10**-23*(x/480*10**6)**7-5.1920577*10**-28*(x/480*10**6)**8
			+(23+273.15))

#Functions to normalize coordinates in video:
def znorm(z):
	return (z*vidz-topz)/(botz-topz)
def xnorm(x):
	return (x-0.5)

#Open tracking file to obtain the length of the video
for file in os.listdir(os.getcwd()):
    if file.endswith("Side_Track.csv"):
        v = file
track = []
with open(v,'r') as csvfile:
	trackreader = csv.reader(csvfile,delimiter = ' ')
	vidtime = len(list(trackreader))/fr
	print "Duration of video: " + str(vidtime) + "s"
	print "Number of frames: " + str(vidtime*fr)
	csvfile.seek(0)
	#In the blender file, X=0 is the center, Z=0 is the top plate.
	i = 0
	for row in trackreader:
		row.pop(0)
		track.append([i/float(fr),znorm(float(row[0])),xnorm(float(row[1]))])
		i = i+1

#Outputs converted tracking data
with open('track_summary.csv', 'wb') as csvfile:
    summwriter = csv.writer(csvfile, delimiter=' ')
    summwriter.writerow(["Time/s","Z","X"])
    for i in track:
		summwriter.writerow(i)

#Converts voltage data to actual values
for file in os.listdir(os.getcwd()):
    if file.endswith(".lvm"):
        v = file
lvm = []
with open(v,'r') as csvfile:
	lvmreader = csv.reader(csvfile, delimiter = '\t')
	n = float(len(list(lvmreader)))
	print "Number of voltage data points: " + str(int(n))
	csvfile.seek(0)
	i = 0
	for row in lvmreader:
		row.pop(0)
		lvm.append([vidtime/n*i,
			pressure(float(row[0])), thermistor(float(row[1])), thermocouple(float(row[2]))])
		i = i+1

#Match pressure data
for file in os.listdir(os.getcwd()):
    if file.endswith("Pressure.csv"):
        v = file
press = []
with open(v, 'r') as csvfile:
	pressreader = csv.reader(csvfile, delimiter='\t')
	n = float(len(list(pressreader)))
	print "Number of pressure data points: " + str(int(n))
	csvfile.seek(0)
	i = 1
	for row in pressreader:
		if i <= 3:
			i = i+1
		else:
			press.append([float(row[0])/1000,float(row[1])])


#Take the average of every n rows in a table.
def avg(table, rate):
	avg = []
	j = 0
	l = len(table[0])
	acc = list(0.0 for i in range(l))
	for row in table:
		if j < rate:
			acc = list(row[i] + acc[i] for i in range(l)) #Keep adding to the accumulator
			j = j+1
		else:
			acc = list(acc[i]/float(rate) for i in range(l)) #Take the average
			avg.append(acc)
			acc = list(0.0 for i in range(l))
			acc = list(row[i] + acc[i] for i in range(l)) #Keep adding to the accumulator
			j = 1
	acc = list(acc[i]/float(j) for i in range(l)) #Take the average of the remaining rows
	avg.append(acc)
	return avg

#Outputs converted voltage data
with open('voltage_summary.csv', 'wb') as csvfile:
    summwriter = csv.writer(csvfile, delimiter=' ')
    summwriter.writerow(["Time/s","Pressure/Torr","Bottom/K","Top/K"])
    for i in lvm:
		summwriter.writerow(i)

#For each tracking data point, find the corresponding temperature and pressure data
def nearest(t,table):
	assert isinstance(t, float)
	minimum = 999999999
	diff = 0.0
	entry = None
	for row in table:
		assert isinstance(row[0],float)
		diff = abs(row[0]-t)
		if diff < minimum:
			entry = row
			minimum = abs(row[0]-t)
	return entry

#Take averages of lvm data
lvm = avg(lvm,10)
new_track = []
for i in range(len(track)):
	if int(i%(len(track)/float(M))) == 0:
		new_track.append(track[i])
track = new_track

#Match track data to lvm data
def match(tb1,tb2):
	matched = []
	i = 0
	for t in tb1:
		matched.append(t + nearest(t[0],tb2))
	return matched

matched = match(match(track,lvm),press)

#Set boundary conditions for plates and outer edges.
def boundaryCondition1(r,z):
    if r<r1*N and z<z1*N:
        T[r][z]=T1
        return 1
    elif r<r2*N and z>z2*N:
        T[r][z]=T2
        return 1
    elif r==N-1:
        T[r][z]=T3
        return 1
    else:
        return 0

def boundaryCondition2():
    for i in range(N):
        T[0][i]=T[1][i]
        T[i][0]=T[i][1]
        T[i][N-1]=T[i][N-2]

def tempweight(h):
	temp_avg = T[0][int(h)]*(int(h)+1-h)+T[0][int(h+1)]*(1-(int(h+1)-h))
	return temp_avg

#Do the simulation!
print
simnum = 1
with open("templist.csv","w") as csvfile:
		templistwriter = csv.writer(csvfile, delimiter=' ')
		templistwriter.writerow(["Time/s","Temperature/K","Gradient/Kmm^-1","Pressure/Torr"])
T=[[150 for x in range(N)] for x in range(N)]
for row in matched:
	print
	print "Simulation "+str(simnum)+"/"+str(M)
	print "Time: " + str(row[0])
	start_time = time.time()
	maxErr=100
	pres = row[8]
	T1 = row[6]
	T2 = row[5]
	print "Top plate temperature: " + str(T1) + "K"
	print "Bottom plate temperature: " + str(T2) + "K"
	print "Pressure (Surrounding): " + str(pres) + " Torr"
	#Apply boundary conditions
	for z in range(N):
		for r in range(N):
			boundaryCondition1(r,z)
	while maxErr>errTop:
		maxErr=0
		for z in range(1,N-1):
			for r in range(1,N-1):
				if boundaryCondition1(r,z)==0:
					Tnew=(T[r-1][z]+T[r+1][z]+T[r][z-1]+T[r][z+1])/4+ \
						 (T[r+1][z]-T[r-1][z])/8/(r-1/2)+ \
						 ((T[r+1][z]-T[r-1][z])**2 +(T[r][z+1]-T[r][z-1])**2)*alpha/16/T[r][z]
					err=abs(Tnew-T[r][z])
					T[r][z]=Tnew
					if err>maxErr:
						maxErr=err
		boundaryCondition2()
	#Find the temperature and gradient at the specified location
	h = (row[1]*(z2-z1)+z1)*N
	x = abs(row[2]*N)
	print "Cell (X): " + str(x)
	print "Cell (Z): " + str(h)
	temp=T[0][int(h)]*(int(h)+1-h)+T[0][int(h+1)]*(1-(int(h+1)-h))
	print "Temperature: "+str(temp)
	tempgrad=((T[0][int(h)+1]-T[0][int(h)])*fac*N+
		(T[0][int(h)]-T[0][int(h-1)])*fac*N)/2
	print "Temp Gradient: " + str(tempgrad)
	pres=pres*math.sqrt(temp/295.0944)
	print "Pressure (Particle location): " + str(pres)
	#Writes result to file
	with open("sim"+str("%03d" % simnum)+".csv", 'wb') as csvfile:
		summwriter = csv.writer(csvfile, delimiter=' ')
		print row[0]
		#Start file with time, grid resolution, gap size, pressure and Z.
		summwriter.writerow([row[0],N,gap,pres,row[1]])
		for r in T:
			summwriter.writerow(r)
	with open("templist.csv","a") as csvfile:
		templistwriter = csv.writer(csvfile, delimiter=' ')
		templistwriter.writerow([row[0],temp, tempgrad,pres])
	simnum = simnum+1
	elapsed_time = time.time() - start_time
	print "Time taken: "+str(elapsed_time)+"s"

elapsed_time = time.time() - init_time
print 
print "Total time taken: "+str(elapsed_time)+"s"