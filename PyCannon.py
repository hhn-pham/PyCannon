import matplotlib.pyplot as plt
import numpy as np
import math
from vpython import *

#EULER-RICHARDSON METHOD ---- 2 Dimensional

def dragModelER(vLaunch,theta,vWindx,vWindy,mass):
    #steps
    N = 12000

    #define delta T and constants
    timeTaken = 1000
    deltaT = timeTaken/N
    g = 10
    uT = (mass*g/0.1)**(1/2)

    #create empty arrays with initial conditions
    Z = np.ones(N)
    X = np.zeros(N)
    Y = np.zeros(N)
    T = np.zeros(N)
    Vz = np.zeros(N)
    Vx = np.zeros(N)
    Vy = np.zeros(N)

    X[0] = 0
    Vx[0] = vLaunch*np.cos(math.radians(theta))
    Y[0] = 0
    Vy[0] = 0
    Z[0] = 0
    Vz[0] = vLaunch*np.sin(math.radians(theta))
    T[0] = 0
    i = 0

    #index arrays
    while (Z[i] >= 0):
        i = i + 1
        T[i]  = i*deltaT
        modVelocity = ((Vx[i-1] - vWindx)**2 + (Vy[i-1] - vWindy)**2 + Vz[i-1]**2)**(1/2)
        ax =  -(g/uT**2)*modVelocity*(Vx[i-1] - vWindx)
        ay =  -(g/uT**2)*modVelocity*(Vy[i-1] - vWindy)                  
        az =  -(g + (g/uT**2)*modVelocity*Vz[i-1])
        vxMid = Vx[i-1] + 0.5*deltaT*ax
        vyMid = Vy[i-1] + 0.5*deltaT*ay
        vzMid = Vz[i-1] + 0.5*deltaT*az
        axMid = -(g/uT**2)*modVelocity*(vxMid - vWindx)
        ayMid = -(g/uT**2)*modVelocity*(vyMid - vWindy)
        azMid = -(g + (g/uT**2)*modVelocity*vzMid)
        X[i] = X[i-1] + vxMid*deltaT
        Y[i] = Y[i-1] + vyMid*deltaT
        Z[i] = Z[i-1] + vzMid*deltaT
        Vx[i] = Vx[i-1] + axMid*deltaT
        Vy[i] = Vy[i-1] + ayMid*deltaT
        Vz[i] = Vz[i-1] + azMid*deltaT
    return T, Vx, Vy, Vz, X, Y, Z, i

T, Vx, Vy, Vz, X, Y, Z, i = dragModelER(240,70,20,20,10)
scene = canvas(autoscale = True, width=X[i] + Y[i] + 100, height=Z[i] + 500, center=vector(0,0,-Z[i]-250), background=color.white, forward=vector(0,1,0), up=vector(0,0,1))
ground = box(pos=vector(0,0,-(Z[i] + 500 + X[i]/1000)), size=vector(X[i] + 100,Y[i]*2 + 100,10), color=color.cyan)
projectile = sphere(canvas = scene, pos=vector(-X[i]/2,0,-(Z[i] + 500)), radius=X[i]/100, color=color.red, make_trail=True)

for j in range(0,i):
    rate(i/T[i])
    projectile.pos.x = X[j] - X[i]/2
    projectile.pos.y = Y[j]
    projectile.pos.z = Z[j] - (Z[i] + 500)

projectile.visible = False
del projectile
    
plt.xlabel(r'$x$ / m')
plt.ylabel(r'$z$ / m')
plt.title('Trajectory in the x-z plane')
plt.plot(X[0:i], Z[0:i], c="coral", linewidth=2)
plt.show()