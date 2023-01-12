import marshal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from functools import partial

G = 6.674 * 1e-11 # N * m^2 / kg^2
M = 1.989 * 1e30 #kg

class Planet():

    def __init__(self, name, mass, coord_init, a, perihelie, e, radius): #a : demi-grand axe, e: excentricité, vitesse en ration par rapport à la vitesse de la Terre
        self.name = name
        
        self.coord = np.array(coord_init)
        
        self.mass = mass
        self.radius = radius
        self.axe = a
        self.perihelie = perihelie
        self.e = e

        self.orbit = []
        self.u = - self.axe + self.perihelie
        self.b = self.axe * np.sqrt(1 - self.e**2)

        self.period = 2 * np.pi * np.sqrt(self.axe**3 / (M + self.mass) / G)
        self.omega = 2 * np.pi / self.period

    def calcOrbit(self, time):

        for i in range(time.size):
            self.coord = np.array([self.u + self.axe * np.cos(self.omega * time[i]), self.b * np.sin(self.omega * time[i])])
            self.orbit.append(self.coord)


    def animate(self, frame, scatter):  

        coord = self.orbit[frame]
        scatter.set_offsets(coord)


mass = {"mercure" : 3.285e23, # kg 
        "venus" : 4.867e24,
        "terre" : 5.972e24,
        "mars" : 6.390e23,
        "jupiter" : 1.898e27,
        "saturne" : 5.683e26,
        "uranus" : 8.681e25,
        "neptune" : 1.024e26}

radius = {"mercure" : 2.439 * 1e6, # m
        "venus" : 6.051 * 1e6,
        "terre" : 6.371 * 1e6,
        "mars" : 3.389 * 1e6,
        "jupiter" :6.991 * 1e7,
        "saturne" : 5.823 * 1e7,
        "uranus" : 2.5362 * 1e7,
        "neptune" : 2.462 * 1e7}

a = {"mercure" : 57.909083e9, # m
    "venus" : 108.210e9,
    "terre" : 149.598e9,
    "mars" : 227.956e9,
    "jupiter" : 778.479e9,
    "saturne" : 1432.041e9,
    "uranus" : 2867.043e9,
    "neptune" : 4514.953e9}

e = {"mercure" : 0.2056, # sans dimension 
        "venus" : 0.0068,
        "terre" : 0.0167,
        "mars" : 0.0935,
        "jupiter" : 0.0487,
        "saturne" : 0.0520,
        "uranus" : 0.0469,
        "neptune" : 0.0097}

perihelie = {"mercure" : 46e9, # m
            "venus" : 107.480e9,
            "terre" : 147.095e9,
            "mars" : 206.650e9,
            "jupiter" : 740.595e9,
            "saturne" : 1357.554e9,
            "uranus" : 2732.696e9,
            "neptune" : 4471.050e9}

coordInit = {"mercure" : [perihelie["mercure"],0],
            "venus" : [perihelie["venus"],0],
            "terre" : [perihelie["terre"],0],
            "mars" : [perihelie["mars"],0],
            "jupiter" : [perihelie["jupiter"],0],
            "saturne" : [perihelie["saturne"],0],
            "uranus" : [perihelie["uranus"],0],
            "neptune" : [perihelie["neptune"],0]}

mercure = Planet("mercure", mass["mercure"], coordInit["mercure"], a["mercure"], perihelie["mercure"], e["mercure"], radius["mercure"])
venus = Planet("venus", mass["venus"], coordInit["venus"], a["venus"], perihelie["venus"], e["venus"], radius["venus"])
terre = Planet("terre", mass["terre"], coordInit["terre"], a["terre"], perihelie["terre"], e["terre"], radius["terre"])
mars = Planet("mars", mass["mars"], coordInit["mars"], a["mars"], perihelie["mars"], e["mars"], radius["mars"])
jupiter = Planet("jupiter", mass["jupiter"], coordInit["jupiter"], a["jupiter"], perihelie["jupiter"], e["jupiter"], radius["jupiter"])
saturne = Planet("saturne", mass["saturne"], coordInit["saturne"], a["saturne"], perihelie["saturne"], e["saturne"], radius["saturne"])
uranus = Planet("uranus", mass["uranus"], coordInit["uranus"], a["uranus"], perihelie["uranus"], e["uranus"], radius["uranus"])
neptune = Planet("neptune", mass["neptune"], coordInit["neptune"], a["neptune"], perihelie["neptune"], e["neptune"], radius["neptune"])

planets = np.array([mercure, venus, terre, mars, jupiter, saturne, uranus, neptune])