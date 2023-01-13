import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from functools import partial
import planetes as pl


G = 6.674*1e-11

figure = plt.figure()

animList = []

T, N = 2.649e9, 10000 # s, n
time = np.linspace(0, T, N)

t = np.linspace(0, 2*np.pi, 100) # x-absc for elliptic trajectory

for planet in pl.planets:
    planet.calcOrbit(time, np.array(planet.coord))
    sc = plt.scatter(planet.coord[0], planet.coord[1])
    animList.append(FuncAnimation(figure, partial(planet.animate, scatter=sc), frames=np.arange(N-1), interval = 17))
    plt.plot(planet.u + planet.axe * np.cos(t) , planet.b * np.sin(t), label=planet.name)

plt.plot(-1.131835012683988e9*1e3, -4.386495519097997e9*1e3, "go")
plt.scatter([0],[0], color="red") # the sun
plt.grid(color='lightgray',linestyle='--')

#plt.legend()
plt.show()