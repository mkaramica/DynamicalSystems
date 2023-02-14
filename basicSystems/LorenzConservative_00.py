#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


rho = 28.0
sigma = 10.0
beta = 8.0 / 3.0



def f(state, t):
    x, y, z = state  # Unpack the state vector
    
    xdot=y
    ydot=-x*z
    zdot=x*y
    
    return xdot, ydot, zdot  # Derivatives

X1=[0,0,0]
X2=[0.3,0,-3.5]
X3=[-0.5,2,0]
X4=[3,2,1]

t = np.arange(0.0, 2.0, 0.0001)

states1 = odeint(f, X1, t)
states2 = odeint(f, X2, t)
states3 = odeint(f, X3, t)
states4 = odeint(f, X4, t)

x1=states1[:,0]
y1=states1[:,1]
z1=states1[:,2]

x2=states2[:,0]
y2=states2[:,1]
z2=states2[:,2]

x3=states3[:,0]
y3=states3[:,1]
z3=states3[:,2]

x4=states4[:,0]
y4=states4[:,1]
z4=states4[:,2]


v1=[]
v2=[]
v3=[]
A=[]

for i in range(len(t)):
    v1.append([x2[i]-x1[i],y2[i]-y1[i],z2[i]-z1[i]])
    v2.append([x3[i]-x1[i],y3[i]-y1[i],z3[i]-z1[i]])
    v3.append([x4[i]-x1[i],y4[i]-y1[i],z4[i]-z1[i]])
    A.append(np.dot(v1[i],np.cross(v2[i],v3[i])))

plt.figure(1)

plt.plot(A)

plt.show()

fig = plt.figure(2)
ax = fig.gca(projection="3d")


ax.plot(x1, y1, z1,'b')
ax.plot(x2, y2, z2,'b')
ax.plot(x3, y3, z3,'b')
ax.plot(x4, y4, z4,'b')

ax.plot(X1[0],X1[1],X1[2],'*r')
ax.plot(X2[0],X2[1],X2[2],'*r')
ax.plot(X3[0],X3[1],X3[2],'*r')
ax.plot(X4[0],X4[1],X4[2],'*r')


ax.plot(x1[-1],y1[-1],z1[-1],'*g')
ax.plot(x2[-1],y2[-1],z2[-1],'*g')
ax.plot(x3[-1],y3[-1],z3[-1],'*g')
ax.plot(x4[-1],y4[-1],z4[-1],'*g')






ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')



plt.draw()
plt.show()