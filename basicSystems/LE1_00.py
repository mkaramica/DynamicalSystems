#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

a=15
b=1

def f(state, t):
    x, y, z = state  # Unpack the state vector
    
    xdot=y
    ydot=-x + y*z
    zdot=-x - a*x*y - b*x*z
    
    return xdot, ydot, zdot  # Derivatives

state0 = [0,0.5,0.5]
t = np.arange(0.0, 400.0, 0.01)

states = odeint(f, state0, t)

x=states[:,0]
y=states[:,1]
z=states[:,2]

fig = plt.figure()
ax = fig.gca(projection="3d")

ax.plot([0,0],[0,0],[min(z),max(z)],'r')
ax.plot(x, y, z,'b')
ax.plot(state0[0],state0[1],state0[2],'*')
#ax.plot(0,0,0,'*r')


ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')



plt.draw()
plt.show()