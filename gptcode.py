import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.674e-11 # gravitational constant

# Initial conditions
m1 = 1.989e30 # mass of first body (e.g. Sun)
m2 = 5.972e24 # mass of second body (e.g. Earth)
r1 = np.array([0, 0]) # initial position of first body
r2 = np.array([149.6e9, 0]) # initial position of second body
v1 = np.array([0, 0]) # initial velocity of first body
v2 = np.array([0, 29800]) # initial velocity of second body

# Time step and number of steps
dt = 3600 # 1 hour
n_steps = 365*24 # 1 year

# Arrays to store positions and velocities
r1_array = np.zeros((n_steps+1, 2))
r2_array = np.zeros((n_steps+1, 2))
v1_array = np.zeros((n_steps+1, 2))
v2_array = np.zeros((n_steps+1, 2))

# Initialize position and velocity arrays
r1_array[0] = r1
r2_array[0] = r2
v1_array[0] = v1
v2_array[0] = v2

# Function to calculate force
def force(r1, r2, m1, m2, G):
    r_vec = r1 - r2
    r_mag = np.linalg.norm(r_vec)
    return -G*m1*m2/r_mag**3*r_vec

# Leapfrog method
for i in range(n_steps):
    # Half step for velocity
    v1_array[i+1] = v1_array[i] + 0.5*dt*force(r1_array[i], r2_array[i], m1, m2, G)
    v2_array[i+1] = v2_array[i] + 0.5*dt*force(r2_array[i], r1_array[i], m2, m1, G)
    # Full step for position
    r1_array[i+1] = r1_array[i] + dt*v1_array[i+1]
    r2_array[i+1] = r2_array[i] + dt*v2_array[i+1]


plt.plot(r1_array[0], r1_array[1])
plt.plot(r2_array[0], r2_array[1])
plt.show()
