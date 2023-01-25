import numpy as np
import matplotlib.pyplot as plt
from voyageSpatial import Planet, Sonde, leapfrog
import planetes as pl

G = 6.674*1e-11

massPlanet = {"Terre" : 1e24, "Lune" : 7.6e24} #kg
dTerreLune = 3.84400 * 1e4 # m
vMax = 1.7 * 1e4 # km/s

day = 86400 # s
hour = 3600 # s
minute = 60 # s
year = 365 * day

terre = Planet("terre", massPlanet["Terre"], [0, 0], pl.terre.radius)

planetArray = np.array([terre])

R = 35000 * 1e3 #m
vyIni = np.sqrt(G * terre.mass / R)
theta = np.pi / 2 # 0 < theta < 2 * np.pi
vDir = np.array([np.cos(theta), -np.sin(theta)])

sonde = Sonde([R, 0], [-1000, vyIni/2], [0, 0])

init = [sonde.x, sonde.y, sonde.vx, sonde.vy]

h = np.array([hour])#, day, 7*day, 14*day, 30*day, year])

#h = np.array([day])
hName = np.array(["hour", "day", "week", "two weeks", "month", "year"])
T = np.arange(1, 15) * year * 3
diffRelatEm = np.empty(T.size)

for k in range(h.size):
    for i in range(T.size):
        print("Caculating for : h =", h[k], "T =", T[i])
        t, u = leapfrog(init, sonde, planetArray, T[i], h[k])

        Em = np.empty(t.size)

        for j in range(u.shape[1]):
            sonde.coord = u[:2, j]
            sonde.speed = u[2:4, j]
            Em[j] = sonde.em(planetArray)

        diffRelatEm[i] = 2 * np.abs(np.max(Em) - np.min(Em)) / np.abs(np.max(Em) + np.min(Em))
    
    plt.plot(T, diffRelatEm, 'o-', label=hName[k])

plt.xlabel("T [year]")
plt.ylabel("Em Difference Relative")
plt.title("dt = " + str(hour/3600) + " heures")
plt.grid()
plt.show()