#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D


mio=3



def f(state, t):
    x, y = state  # Unpack the state vector
    
    xdot=y
    ydot=mio*y*(1-x**2)-x 
        
    
    return xdot, ydot  # Derivatives

def f_rotated(state, t):    # 90 degree counter-clockwise
    x, y = state  # Unpack the state vector
    
    xdot= mio*x*(1-y**2)+y  
    ydot=-x 
    
    return xdot, ydot  # Derivatives   
    
x0=.2
y0=.2


state0= [x0,y0]


t = np.arange(0.0, 200.0, 0.01)

states = odeint(f, state0, t)

x=states[:,0]
y=states[:,1]

states_rotated = odeint(f_rotated, state0, t)

x_rotated=states_rotated[:,0]
y_rotated=states_rotated[:,1]


fig = plt.figure(1)
ax = fig.gca()
ax.plot(x,y,'r')
ax.plot(x0,y0,'*r')
ax.plot(x_rotated,y_rotated,'.b')
ax.plot(y,-x,'g')

plt.gca().set_aspect('equal', adjustable='box')

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.grid(True)
plt.draw()
plt.show()