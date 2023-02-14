#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

n0AVG=10000

def f(state, t):
    x, y, z = state  # Unpack the state vector
    
    xdot = y
    ydot = -x +  y * z
    zdot = 1 - y**2
    
    return xdot, ydot, zdot  # Derivatives



x0Mat=[0.1,-1.1,-1,0,.2,.3,1,2,0,-1,.2,-.3]
y0Mat=[0.1,-1.1,3,2,.3,-0.2,0,.2,-0.5,-1.5,-.8,-.3]
z0Mat=[0.1,-1.1,1,0.1,.2,-.3,-1,0,0.2,2,0,-.25]

t = np.arange(0.0, 400.0, 0.01)


Len_t=len(t)

print("Len_t=",Len_t)
print("---------------")

'''
fig = plt.figure(1)
ax = fig.gca(projection="3d")
'''

fig = plt.figure(1)
ax = fig.add_subplot(2, 2, 1, projection='3d')

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')


indxSubPlot=1

zAVGmat=[]


for i in range(0,12):
    state0 = [x0Mat[i], y0Mat[i], z0Mat[i]]
    states = odeint(f, state0, t)  
    
    zVec=states[n0AVG-1:-1,2]
    zAVG=sum(zVec)/(Len_t-n0AVG)
    zAVGmat.append(zAVG)
    
    print(i,zAVG)
    
    
    
    ax.plot(states[1500:-1, 0], states[1500:-1, 1], states[1500:-1, 2],label=str(i))
    if i==2 or i==5 or i==8:
        ax.legend(loc = 'best')
        indxSubPlot=indxSubPlot+1
        ax = fig.add_subplot(2, 2, indxSubPlot, projection='3d') 
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

print("---------------")
print("Total Average of z:" ,sum(zAVGmat)/12)        
ax.legend(loc = 'best')
plt.draw()
plt.show()


