import numpy as np
import matplotlib.pyplot as plt

M = 5.972e24
Ms = 1.989 * 1e30 # kg
TS = 149.598e9 # m
m = 100
G = 6.674*1e-11

minute = 60 # s
hour = 60 * minute
day = 24 * hour
week = 7 * day
month = 30.5 * day
year = 365 * day

def calculateEm(r, v, rPlanet):
    return 1/2 * m * np.linalg.norm(v)**2 - G * M * m / np.linalg.norm(rPlanet - r)

def a(sondeCoord, planetCoord):
    return -G * M * (planetCoord - sondeCoord)/ np.linalg.norm(planetCoord - sondeCoord)**3

def leapfrog(u_ini, planetArray, t, h):
    N = t.size
    # Initialisation du tableau
    u = np.empty((4, 1)) # Condition initiale
    u[:, 0] = u_ini    

    # Coefficient d'intégration
    w0 = - 2**(1/3) / (2 - 2**(1/3))
    w1 = 1 / (2 - 2**(1/3))

    c = np.array([w1 / 2, (w0 + w1) / 2, (w0 + w1) / 2, w1/2])
    d = np.array([w1, w0, w1])

    em = np.array([0])
    for planet in planetArray:
        em[0] += calculateEm(u[:2, 0], u[2:4, 0], planet.orbit[0, :])

    time = np.array([0])
    i = 0
    T = t[-1]
    timeTracker = 0
    while timeTracker <= T:

        h = np.floor(variableTimestep(i, h, u[:2, i], planetArray, month, minute, 1.1, 0.7))
        time = np.append(time, h + time[-1])
        #4th order Yoshida integrator
        
        r1 = u[:2, i] + c[0] * u[2:4, i] * h

        a1 = 0
        for planet in planetArray:
            a1 += a(r1, planet.orbit[i, :])

        v1 = u[2:4, i] - d[0] * a1 * h

        r2 = r1 + c[1] * v1 * h

        a2 = 0
        for planet in planetArray:
            a2 += a(r2, planet.orbit[i, :])

        v2 = v1 - d[1] * a2 * h

        r3 = r2 + c[2] * v2 * h

        a3 = 0
        for planet in planetArray:
            a3 += a(r3, planet.orbit[i, :])

        v3 = v2 - d[2] * a3 * h

        # print("r3,v3",r3,v3)
        r4 = r3 + c[3] * v3 * h
        v4 = v3 

        x = r4[0]
        y = r4[1]
        vx = v4[0]
        vy = v4[1]

        uNew = np.array([[x], [y], [vx], [vy]])
        u = np.append(u, uNew, axis=1)

        emTemp = 0
        for planet in planetArray:
            # print(planet.orbit[i+1, :], u[:2, i+1])
            emTemp += calculateEm(u[:2, i+1], u[2:4, i+1], planet.orbit[i+1, :])
        em = np.append(em, emTemp)
    
        print(timeTracker, "/", T, end="\r")

        if np.linalg.norm(planetCoord - u[:2, i+1]) <= MAX_RADIUS:
            print("\nError : collision with planet's surface")
            u[:2, i+1] = u[:2, i]
            return t, u, em

        timeTracker += h
        i += 1

    print("\n")
    return time, u, em

def variableTimestep(i, h, sondeCoord, planetArray, hMax, hMin, upScale, downScale, hillSphereScaleFactor=3):
    newh = h
    for planet in planetArray:
        threshold = hillSphereScaleFactor * planet.hillSphereRadius

        if np.linalg.norm(planet.orbit[i, :] - sondeCoord) <= threshold:
            newh *= downScale
            if newh < hMin:
                newh = hMin
        else:
            newh *= upScale
            if newh > hMax:
                newh = hMax
    return newh
            
MAX_RADIUS = 6.3e6

class Planet():

    def __init__ (self, mass, orbit, hillSphereRadius, radius=0):
        self.mass=mass
        self.orbit=orbit
        self.raidus=radius
        self.hillSphereRadius = hillSphereRadius

rayonSphereHill = TS * (1-0.0167) * (M / Ms / 3)**(1/3)# m
planetCoord = np.array([0,0])
planet2Coord = np.array([-2.603e11, 3.9e11])

vitesseSonde = 1666 # m/s
sondeInitCoord = np.array([-5e8, -1e9]) #m
sondeInitSpeed = np.array([1, 2]) / np.linalg.norm(np.array([1, 2])) * vitesseSonde

vitesseTerre = 2978 # m/s
timeImpact = np.linalg.norm(sondeInitCoord) / vitesseSonde
distanceInitTerre = timeImpact * vitesseTerre
u_ini = np.append(sondeInitCoord, sondeInitSpeed)

T = month
dt = 1 #10 * minute

t = np.arange(0, T, dt)

P = year # period 
R = 2e9
# planet_orbit = np.array([np.cos(2 * np.pi * t / P), np.sin(2 * np.pi * t / P)]).T * R
# print("orbit\n", planet_orbit)
# planet_orbit = np.zeros([t.size, 2])
planet_orbit = np.array([np.linspace(-distanceInitTerre, t[-1]*vitesseTerre, t.size), np.zeros(t.size)]).T
# print(planet_orbit)
planet_orbit2 = np.ones([t.size, 2]) * planet2Coord

planet = Planet(M, planet_orbit, rayonSphereHill)
planet2 = Planet(M, planet_orbit2, rayonSphereHill)

planetArray = np.array([planet,])
time,u,em = leapfrog(u_ini, planetArray, t, dt)


# Sonde
plt.plot(u[0, :], u[1, :], '-', label = "Trajectoire de la sonde")
plt.plot(u[0, 0], u[1, 0], 'bo', label = "Position initiale de la sonde")
# plt.plot(u[0, -1], u[1, -1], 'b*', label = "sonde end")

# Planet
for planet in planetArray:
    plt.plot(planet_orbit[::100, 0], planet_orbit[::100, 1], label = "Trajectoire de la Terre")
    plt.plot(planet.orbit[0, 0], planet.orbit[0, 1], 'ro', label = "Terre")
    # plt.plot(MAX_RADIUS * np.cos(np.linspace(0, 2 * np.pi, 1000)) + planet.orbit[0, 0], MAX_RADIUS * np.sin(np.linspace(0, 2 * np.pi, 1000)) + planet.orbit[0, 1], 'r--', label="rayon")

    plt.plot(planet_orbit[-1, 0], planet_orbit[-1, 1], 'r*', label = "Position finale de la planète")

# # Planet2
# # plt.plot(planet_orbit2[:, 0], planet_orbit2[:, 1], label = "planet2")
# plt.plot(planet_orbit2[0, 0, 0], planet_orbit2[0, 0, 1], 'go', label = "planet2 start")
# # plt.plot(planet_orbit2[-1, 0], planet_orbit2[-1, 0], 'r*', label = "planet2 end")

# plt.plot(planetCoord[0], planetCoord[1], 'ro', label='planet')
# plt.plot(rayonSphereHill * np.cos(np.linspace(0, 2 * np.pi, 1000)) + planetCoord[0], rayonSphereHill * np.sin(np.linspace(0, 2 * np.pi, 1000)) + planetCoord[1], 'r--')
# plt.plot(rayonSphereInfluence * np.cos(np.linspace(0, 2 * np.pi, 1000)) + planet2Coord[0], rayonSphereInfluence * np.sin(np.linspace(0, 2 * np.pi, 1000)) + planet2Coord[1], 'g--')

plt.xlim(-1.2e9, 1.2e9)
plt.ylim(-1.2e9, 1.2e9)
plt.xlabel("x [m]")
plt.ylabel("y [m]")

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.legend()
plt.grid()
plt.title("Fronde gravitationnelle dans le référenciel planétaire")
plt.show()

# Em
emRelatif = np.abs(2 * (np.max(em) - np.min(em)) / (np.max(em) + np.min(em)))
# plt.plot(time, em, '-', label="em")
print("Em =", emRelatif)

# Sonde speed
print(u[2:4], u[2, :], u[3, :])
vx = u[2, :]
vy = u[3, :]
v = np.sqrt(np.power(vx, 2) + np.power(vy, 2))
print(v)

plt.plot(time, v, '-')
plt.xlabel("t [s]")
plt.ylabel("v [m/s]")

plt.title("Évolution temporelle de la vitesse sur un an")
plt.grid()
plt.legend()
plt.show()