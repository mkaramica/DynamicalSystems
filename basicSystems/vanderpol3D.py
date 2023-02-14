#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D


mio=3



def f_2D(state, t):
    x, y = state  # Unpack the state vector
    
    xdot=y
    ydot=mio*y*(1-x**2)-x    
    
    return xdot, ydot  # Derivatives

def f_3D(state, t):
    x, y, z = state  # Unpack the state vector
    
    xdot=y
    ydot=z-x    
    zdot=mio*(z-x)*(1-x**2)-2*mio*x*y**2
    
    return xdot, ydot, zdot  # Derivatives
    
    
x0=2
y0=2
z0=mio*y0*(1-x0**2)

state0_2d = [x0,y0]
state0_3d = [x0,y0,z0]

t = np.arange(0.0, 200.0, 0.01)

states_2d = odeint(f_2D, state0_2d, t)

x_2d=states_2d[:,0]
y_2d=states_2d[:,1]

z_2d=[mio*y_2d[i]*(1-x_2d[i]**2)-x_2d[i] for i in range(len(x_2d))]

states_3d = odeint(f_3D, state0_3d, t)

x_3d=states_3d[:,0]
y_3d=states_3d[:,1]
z_3d=states_3d[:,2]



fig = plt.figure(1)
ax = fig.gca(projection="3d")
ax.plot(x_2d,y_2d,z_2d,'r')
ax.plot(x0,y0,z0,'*r')
ax.plot(x_3d,y_3d,z_3d,'b')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
ax.grid(True)
ax.view_init(-130, 130)
plt.draw()
plt.show()