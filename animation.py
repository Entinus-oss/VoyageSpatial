import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from functools import partial
import planetes as pl


G = 6.674*1e-11 # constante gravitaionnel

figure = plt.figure()

animList = []

T, N = 2.649e9, 10000 # s, n
time = np.linspace(0, T, N)

t = np.linspace(0, 2*np.pi, 100) # x-absc for elliptic trajectory

### Calcule la trajectoire et anime l'orbite des planètes ###
for planet in pl.planets:
    planet.calcOrbit(time)
    sc = plt.scatter(planet.coord[0], planet.coord[1]) # plot la position intiale planètes
    animList.append(FuncAnimation(figure, partial(planet.animate, scatter=sc), frames=np.arange(N-1), interval = 30)) # anime et rafraichi toute les 30 ms
    plt.plot(planet.u + planet.axe * np.cos(t) , planet.b * np.sin(t), label=planet.name) # plot l'orbite des planètes

plt.scatter([0],[0], color="red") # affiche le soleil
plt.grid(color='lightgray',linestyle='--') 

#plt.legend()
plt.show()

