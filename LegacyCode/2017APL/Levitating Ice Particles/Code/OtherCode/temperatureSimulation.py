import numpy as np

fac = 0.6/12.7	#Scaling factor
r1 = 25.4/2		#Top plate radius
r2 = 33.7/2		#Bottom plate radius
gap = 10.36		#Gap size
N = 128		#Grid resolution of simulation
M = 100			#Number of simulations
T1 = 77			#Temp of top plate
T2 = 277		#Temp of bottom plate
T3 = 298		#Temp of surroundings
errTop = 0.01	#Maximum error; the smaller the slower the simulation
r1 = r1*fac		#N2 plate radius
z1 = 0.094		#Top plate
r2 = r2*fac		#Cu plate radius
z2 = z1+gap*fac	#Bottom plate height
err = 100
alpha = 0.918

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
        
T=[[150 for x in range(N)] for x in range(N)]
maxErr = 300
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
np.savetxt("C:/Users/Admin/Desktop/PaperData/Data/temperaturedata.txt",T,fmt='%.2f')