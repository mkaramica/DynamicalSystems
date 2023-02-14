#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

# Paper: On a 3-D generalized Hamiltonian model with conservative and dissipative chaotic flows

def f(state, t):
    
    a=1
    b=0
    c=0
    u=1
    
    
    
    x, y, z = state  # Unpack the state vector
    
    xdot = -c * x + a * y
    ydot = -a * x -  y * z
    zdot = -b * z +y ** 2 - u
    
    return xdot, ydot, zdot  # Derivatives


fig = plt.figure()
ax = fig.gca(projection="3d")


state0 = [0, -2, 0]
t = np.arange(0.0, 200.0, 0.01)

states = odeint(f, state0, t)


ax.plot(states[1000:-1, 0], states[1000:-1, 1], states[1000:-1, 2],'r')
ax.plot(state0[0],state0[1],state0[2],'*r')


state0 = [2, -2, -2]
states = odeint(f, state0, t)
t = np.arange(0.0, 200.0, 0.01)
ax.plot(states[1000:-1, 0], states[1000:-1, 1], states[1000:-1, 2],'b')
ax.plot(state0[0],state0[1],state0[2],'*b')



ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')


plt.draw()
plt.show()