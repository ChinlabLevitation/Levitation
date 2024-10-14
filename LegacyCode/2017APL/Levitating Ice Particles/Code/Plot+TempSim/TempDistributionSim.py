import math

# N2 diameter = 25.4mm, Cu diameter = 33.8 mm, separation = 10.5 mm
# all coordinates scaled by 0.6/12.7 mm.  

N=128
T1=77
T2=477
T3=290
r1=0.6
z1=0.251
r2=0.796
z2=0.748
errTop=0.01
err=100
alpha=0.918
maxErr=100
T=[[150 for x in range(N)] for x in range(N)]
dataOutput=open("C:/Users/Admin/Desktop/Levitating Ice Particles/Code/Precise.txt","w")

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

for z in range(N):
    for r in range(N):
        boundaryCondition1(r,z)

while maxErr>errTop:
    maxErr=0
    for z in range(1,N-1):
        for r in range(1,N-1):
            if boundaryCondition1(r,z)==0:
                Tnew=1.0*(T[r-1][z]+T[r+1][z]+T[r][z-1]+T[r][z+1])/4+ \
                     (T[r+1][z]-T[r-1][z])/8/(r-1/2)+ \
                     ((T[r+1][z]-T[r-1][z])**2 +(T[r][z+1]-T[r][z-1])**2)*alpha/16/T[r][z]
                err=abs(Tnew-T[r][z])
                T[r][z]=Tnew
                if err>maxErr:
                    maxErr=err
    boundaryCondition2()

for z in range(N):
    for r in range(N):
        #print(T[r][z]," ",end="", file=dataOutput)
        dataOutput.write('%f ' % T[r][z])
#dataOutput.write('%d' % maxErr)

dataOutput.close()
