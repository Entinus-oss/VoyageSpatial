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

def find_intersection(trajectory1, trajectory2, radius):
    intersections = []
    for point1 in trajectory1:
        for point2 in trajectory2:
            distance = np.linalg.norm(point1 - point2)
            if distance <= radius:
                intersections.append((point1, point2))
    return intersections

MAX_RADIUS = 6.3e5 # m

T = 1 * year
dt = 10 * minute

t = np.arange(0, T, dt)
P = year * 3
TS = 150e9
R = 35000 * 1e3 # m
vy0 = np.sqrt(G * M / R)

planet_orbit = np.array([TS * np.cos(2 * np.pi * t / P), TS * np.sin(2 * np.pi * t / P)])
# planet_orbit = np.array([np.linspace(0, 1e6, t.size), np.linspace(0, 1e6, t.size)])

u_ini =[R + TS, 0, 0, vy0 + TS * 2 * np.pi / P]
t, u, em = leapfrog(u_ini, planet_orbit, t, dt)

emRelatif = np.abs(2 * (np.max(em) - np.min(em)) / (np.max(em) + np.min(em)))

# Sonde
plt.plot(u[0, :], u[1, :], label = "sonde")
plt.plot(u[0, 0], u[1, 0], 'bo', label = "sonde start")
plt.plot(u[0, -1], u[1, -1], 'b*', label = "sonde end")

# Planet
plt.plot(planet_orbit[0, :], planet_orbit[1, :], label = "planet")
plt.plot(planet_orbit[0, 0], planet_orbit[1, 0], 'ro', label = "planet start")
plt.plot(planet_orbit[0, -1], planet_orbit[1, -1], 'r*', label = "planet end")

plt.legend()
plt.show()

# Em
plt.plot(t, em, '-', label="em")
print(emRelatif)

plt.legend()
plt.show()
