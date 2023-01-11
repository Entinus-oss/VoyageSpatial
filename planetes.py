import marshal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from functools import partial


class Planet():

    def __init__(self, name, m, coord_init, a=0, perihelie=0, e=0): #a : demi-grand axe, e: excentricité, vitesse en ration par rapport à la vitesse de la Terre
        self.m = m
        self.name = name
        #self.radius = radius
        self.coord = np.array(coord_init)
        self.a = a
        self.perihelie = perihelie
        #self.aphelie = aphelie
        self.e = e
        #self.v = v
        self.u = -self.a+ self.perihelie
        self.b = self.a*np.sqrt(1-self.e**2)

    def actualisation_position(self,t, scatter):
        self.coord = [self.u+self.a*np.cos(t), self.b*np.sin(t)]
        scatter.set_offsets(self.coord)


M = 5.972e24 
mass = {"mercure" : 3.285e23,
        "venus" : 4.867e24,
        "terre" : 5.972e24,
        "mars" : 6.390e23,
        "jupiter" : 1.898e27,
        "saturne" : 5.683e26,
        "uranus" : 8.681e25,
        "neptune" : 1.024e26,   
        "soleil" : 1.989e30}
a = {"mercure" : 57.909083e6,
    "venus" : 108.210e6,
    "terre" : 149.598e6,
    "mars" : 227.956e6,
    "jupiter" : 778.479e6,
    "saturne" : 1432.041e6,
    "uranus" : 2867.043e6,
    "neptune" : 4514.953e6}
e = {"mercure" : 0.2056,
        "venus" : 0.0068,
        "terre" : 0.0167,
        "mars" : 0.0935,
        "jupiter" : 0.0487,
        "saturne" : 0.0520,
        "uranus" : 0.0469,
        "neptune" : 0.0097}
perihelie = {"mercure" : 46e6,
            "venus" : 107.480e6,
            "terre" : 147.095e6,
            "mars" : 206.650e6,
            "jupiter" : 740.595e6,
            "saturne" : 1357.554e6,
            "uranus" : 2732.696e6,
            "neptune" : 4471.050e6}
coordInit = {"mercure" : [perihelie["mercure"],0],
            "venus" : [perihelie["venus"],0],
            "terre" : [perihelie["terre"],0],
            "mars" : [perihelie["mars"],0],
            "jupiter" : [perihelie["jupiter"],0],
            "saturne" : [perihelie["saturne"],0],
            "uranus" : [perihelie["uranus"],0],
            "neptune" : [perihelie["neptune"],0]}

mercure = Planet("mercure", mass["mercure"], coordInit["mercure"], a["mercure"], perihelie["mercure"], e["mercure"])
venus = Planet("venus", mass["venus"], coordInit["venus"], a["venus"], perihelie["venus"], e["venus"])
terre = Planet("terre", mass["terre"], coordInit["terre"], a["terre"], perihelie["terre"], e["terre"])
mars = Planet("mars", mass["mars"], coordInit["mars"], a["mars"], perihelie["mars"], e["mars"])
jupiter = Planet("jupiter", mass["jupiter"], coordInit["jupiter"], a["jupiter"], perihelie["jupiter"], e["jupiter"])
saturne = Planet("saturne", mass["saturne"], coordInit["saturne"], a["saturne"], perihelie["saturne"], e["saturne"])
uranus = Planet("uranus", mass["uranus"], coordInit["uranus"], a["uranus"], perihelie["uranus"], e["uranus"])
neptune = Planet("neptune", mass["neptune"], coordInit["neptune"], a["neptune"], perihelie["neptune"], e["neptune"])