# orbit around rectilinear motion
length = 1e9 # m 
planet_orbit = np.array([np.linspace(0, length, t.size), np.linspace(0, length, t.size)])

# orbit around stationnary
P = year # period 
planet_orbit = np.array([np.cos(2 * np.pi * t / P), np.sin(2 * np.pi * t / P)])