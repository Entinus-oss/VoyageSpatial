import numpy as np
import matplotlib.pyplot as plt

class Planet():

    def __init__(self, mass, coord, orbit, radius):
        self.mass = mass
        self.coord = np.array(coord)
        self.orbit = orbit
        self.radius = radius

class Sonde():

    def __init__(self, mass, coord, speed):
        self.mass = mass
        self.coord = np.array(coord)
        self.speed = np.array(speed)
        self.orbit = []

def a(sondeIntermediateCoord, PlanetArray):
    acc = 0
    for planet in PlanetArray:
        acc += -G * planet.mass * (planet.coord - sondeIntermediateCoord)/ np.linalg.norm(planet.coord - sondeIntermediateCoord)**3
    return acc

def leapfrog(Sonde, PlanetArray, dt):

    # Coefficient d'intégration
    w0 = - 2**(1/3) / (2 - 2**(1/3))
    w1 = 1 / (2 - 2**(1/3))

    c = np.array([w1 / 2, (w0 + w1) / 2, (w0 + w1) / 2, w1/2])
    d = np.array([w1, w0, w1])

    #4th order Yoshida integrator
    
    r1 = Sonde.coord + c[0] * Sonde.coord * dt
    v1 = Sonde.speed + d[0] * a(r1, PlanetArray) * dt

    r2 = r1 + c[1] * v1 * dt
    v2 = v1 + d[1] * a(r2, PlanetArray) * dt

    r3 = r2 + c[2] * v2 * dt
    v3 = v2 + d[2] * a(r3, PlanetArray) * dt

    Sonde.coord = r3 + c[3] * v3 * dt
    Sonde.speed = v3

    Sonde.orbit.append(Sonde.coord)


G = 6.674*1e-11

Mt = 5.972e24
Ms = 1.989e30

Rt = 6e5
Rs = 7e8

m = 100
RStationnaryOrbit = 3.5e7

minute = 60 # s
hour = 60 * minute
day = 24 * hour
week = 7 * day
year = 365 * day


sonde = Sonde(m, [RStationnaryOrbit, 0], [0, 0])
terre = Planet(Ms, [0, 0], [], 1e4)

planetArray = np.array([terre])

T = 1 * year
dt = 1/3 * hour

t = np.arange(0, T, dt)

N = t.size

for i in range(N):
    leapfrog(sonde, planetArray, dt)

plt.plot(terre.coord[0], terre.coord[1], 'ro')
plt.plot(t, sonde.orbit)
plt.show()
