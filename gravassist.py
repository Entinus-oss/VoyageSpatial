import numpy as np
import matplotlib.pyplot as plt

M = 5.972e24
m = 100
G = 6.674*1e-11

minute = 60 # s
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365 * day

def calculateEm(r, v, rPlanet):
    return 1/2 * m * np.linalg.norm(v)**2 - G * M * m / np.linalg.norm(rPlanet - r)

def a(sondeCoord, planetCoord):
    return -G * M * (planetCoord - sondeCoord)/ np.linalg.norm(planetCoord - sondeCoord)**3

def leapfrog(u_ini, planet_orbit, t, h):
    N = t.size
    # Initialisation du tableau
    u = np.empty((4, N)) # Condition initiale
    u[:, 0] = u_ini    

    # Coefficient d'intégration
    w0 = - 2**(1/3) / (2 - 2**(1/3))
    w1 = 1 / (2 - 2**(1/3))

    c = np.array([w1 / 2, (w0 + w1) / 2, (w0 + w1) / 2, w1/2])
    d = np.array([w1, w0, w1])

    em = np.empty(t.size)
    em[0] = calculateEm(u[:2, 0], u[2:4, 0], planet_orbit[:, 0])

    for i in range(N-1):

        planetCoord = planet_orbit[:, i]

        #4th order Yoshida integrator
        
        r1 = u[:2, i] + c[0] * u[2:4, i] * h
        v1 = u[2:4, i] - d[0] * a(r1, planetCoord) * h

        r2 = r1 + c[1] * v1 * h
        v2 = v1 - d[1] * a(r2, planetCoord) * h

        r3 = r2 + c[2] * v2 * h
        v3 = v2 - d[2] * a(r3, planetCoord) * h

        u[:2, i+1] = r3 + c[3] * v3 * h
        u[2:4, i+1] = v3 

        em[i+1] = calculateEm(u[:2, i+1], u[2:4, i+1], planet_orbit[:, i+1])
    
        print(i, "/", N-2, end="\r")

        if np.linalg.norm(planetCoord - u[:2, i+1]) <= MAX_RADIUS:
            print("Error : collision with planet's surface")
            u[:2, i+1] = u[:2, i]
            return t, u, em

    return t, u, em

MAX_RADIUS = 6.3e5

rayonSphereInfluence = 1e9 # m
planetCoord = np.array([0,0])
sondeInitCoord = np.array([0.5e9, -2e9]) #m
sondeInitSpeed = np.array([0, 6e3 / 3600 * 1e3]) # m/s

u_ini = np.append(sondeInitCoord, sondeInitSpeed)

T = 100*year
dt = hour

t = np.arange(0, T, dt)

planet_orbit = np.array([np.zeros([t.size, 2])])
t,u,em = leapfrog(u_ini, planet_orbit, t, dt)

plt.plot(planetCoord[0], planetCoord[1], 'ro', label='planet')
plt.plot(rayonSphereInfluence * np.cos(np.linspace(0, 2 * np.pi, 1000)) - planetCoord[0], rayonSphereInfluence * np.sin(np.linspace(0, 2 * np.pi, 1000)) - planetCoord[1], '--')

# Sonde
plt.plot(u[0, :], u[1, :], label = "sonde")
plt.plot(u[0, 0], u[1, 0], 'bo', label = "sonde start")
plt.plot(u[0, -1], u[1, -1], 'b*', label = "sonde end")

# Planet
plt.plot(planet_orbit[:, 0], planet_orbit[:, 1], label = "planet")
plt.plot(planet_orbit[0, 0], planet_orbit[0, 1], 'ro', label = "planet start")
plt.plot(planet_orbit[-1, 0], planet_orbit[-1, 0], 'r*', label = "planet end")

plt.legend()
plt.show()

# Em
emRelatif = np.abs(2 * (np.max(em) - np.min(em)) / (np.max(em) + np.min(em)))
plt.plot(t, em, '-', label="em")
print(emRelatif)

plt.legend()
plt.show()
