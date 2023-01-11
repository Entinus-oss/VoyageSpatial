import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from functools import partial
import planetes as pl


G = 6.674*1e-11

figure = plt.figure()
sc = plt.scatter(pl.mercure.coord[0], pl.mercure.coord[1])
sc2 = plt.scatter(pl.venus.coord[0], pl.venus.coord[1])
sc3 = plt.scatter(pl.terre.coord[0], pl.terre.coord[1])
sc4 = plt.scatter(pl.mars.coord[0], pl.mars.coord[1])
sc5 = plt.scatter(pl.jupiter.coord[0], pl.jupiter.coord[1])
sc6 = plt.scatter(pl.uranus.coord[0], pl.uranus.coord[1])
sc7 = plt.scatter(pl.neptune.coord[0], pl.neptune.coord[1])

animation = FuncAnimation(figure,partial(pl.mercure.actualisation_position, scatter = sc),frames = np.arange(0, 2*np.pi, 0.1),interval = 10)
animation2 = FuncAnimation(figure, partial(pl.venus.actualisation_position, scatter = sc2), frames = np.arange(0, 2*np.pi, 0.1),interval = 10 )
animation3 = FuncAnimation(figure,partial(pl.terre.actualisation_position, scatter = sc3),frames = np.arange(0, 2*np.pi, 0.1),interval = 10)
animation4 = FuncAnimation(figure, partial(pl.mars.actualisation_position, scatter = sc4), frames = np.arange(0, 2*np.pi, 0.1),interval = 10 )
animation5 = FuncAnimation(figure,partial(pl.jupiter.actualisation_position, scatter = sc5),frames = np.arange(0, 2*np.pi, 0.1),interval = 10)
animation6 = FuncAnimation(figure, partial(pl.uranus.actualisation_position, scatter = sc6), frames = np.arange(0, 2*np.pi, 0.1),interval = 10 )
animation7 = FuncAnimation(figure, partial(pl.neptune.actualisation_position, scatter = sc7), frames = np.arange(0, 2*np.pi, 0.1),interval = 10 )

t = np.linspace(0, 2*np.pi, 100)
fig = plt.figure()
ax = fig.add_subplot()
plt.plot(pl.mercure.u+pl.mercure.a*np.cos(t) , pl.mercure.b*np.sin(t), label="mercure")
plt.plot(pl.venus.u+pl.venus.a*np.cos(t) , pl.venus.b*np.sin(t), label="venus")
plt.plot(pl.terre.u+pl.terre.a*np.cos(t) , pl.terre.b*np.sin(t), label="terre")
plt.plot(pl.mars.u+pl.mars.a*np.cos(t) , pl.mars.b*np.sin(t), label="mars")
plt.plot(pl.jupiter.u+pl.jupiter.a*np.cos(t) , pl.jupiter.b*np.sin(t), label="jupiter")
plt.plot(pl.uranus.u+pl.uranus.a*np.cos(t) , pl.uranus.b*np.sin(t), label="uranus")
plt.plot(pl.neptune.u+pl.neptune.a*np.cos(t) , pl.neptune.b*np.sin(t), label="neptune")

plt.scatter([0],[0], color="red")
plt.grid(color='lightgray',linestyle='--')

ax.set_aspect('equal')
plt.legend()
plt.show()

