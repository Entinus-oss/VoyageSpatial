import numpy as np
import matplotlib.pyplot as plt
import planetes as pl

G = 6.6e-11 * 1e-3 ** 2 # N * m^2 / kg^2
r = np.linspace(100, 10e10, int(10e6)) # rayons tests
m = 1000 # kg masse de la sonde

planetArray = np.array([pl.mercure, pl.venus, pl.terre, pl.mars, pl.jupiter, pl.saturne, pl.uranus, pl.neptune])
def gravForces(planet):
    forces =  planet.m * G * m / r**2
    plt.plot(r, forces, label = planet.name)

for i in range(planetArray.size):
    gravForces(planetArray[i])

plt.legend()
plt.show()