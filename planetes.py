import marshal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from functools import partial

G = 6.674*1e-11

class Planet():

    def __init__(self, m, coord_init, a=0, perihelie=0, e=0): #a : demi-grand axe, e: excentricité, vitesse en ration par rapport à la vitesse de la Terre
        self.m = m
        self.coord = np.array(coord_init)
        self.axe = a
        self.perihelie = perihelie
        self.e = e
        self.orbit = []

        self.u = -self.axe + self.perihelie
        self.b = self.axe * np.sqrt(1 - self.e**2)
        self.period = 2 * np.pi * np.sqrt(self.axe**3 / self.m / G)
        self.omega = 2 * np.pi / self.period

    # def planetVelocity(self, time): # time un array des temps 
    #     vitesse = np.ones(2)
        
    #     self.coord = [self.u + self.axe * np.cos(time), self.b * np.sin(time)]
    #     radius = np.linalg.norm(self.coord)
    #     velocityScalar = np.sqrt(G * self.m * ((2 * self.axe - radius)/(self.axe * radius)))

    #     vitesse[0] = velocityScalar * -self.axe * np.sin(time)
    #     vitesse[1] = velocityScalar * np.cos(time)
    #     vitesse = vitesse/np.linalg.norm(vitesse)

    #     return vitesse

    # def actualisation_position0(self, time, scatter):
    #     # Calcule les coordonnées grâce à la vitesse 
         
    #     actualisationVitesse = self.planetVelocity(time)

    #     self.coord = [self.coord[0] + actualisationVitesse[0], self.coord[1] + actualisationVitesse[1]]
    #     scatter.set_offsets(self.coord)

    def actualisation_position(self, time, T, N):
        # Calcule les coordonnées grâce à la vitesse 
        #self.perimeter = np.pi*(3*(self.axe+self.b) - np.sqrt(3*self.axe + self.b)*(self.axe * 3*self.b))
        time = np.linspace(0, T, N)
        for i in range(time.size):
            print("time", time[i])
            self.coord = np.array([self.u + self.axe * np.cos(self.omega * time[i]), self.b * np.sin(self.omega * time[i])])
            self.orbit.append(self.coord)
            print("coord", self.coord)


    def animate(self, frame, scatter):  
        # print(frame)
        # print(len(self.orbit))
        coord = self.orbit[frame]
        #print("coord", coord)
        scatter.set_offsets(coord)
        #return scatter


M = 5.972e24 
mass = {"mercure" : 3.285e23, # kg 
        "venus" : 4.867e24,
        "terre" : 5.972e24,
        "mars" : 6.390e23,
        "jupiter" : 1.898e27,
        "saturne" : 5.683e26,
        "uranus" : 8.681e25,
        "neptune" : 1.024e26}
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
perihelie = {"mercure" : 46e9, # km
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

mercure = Planet(mass["mercure"], coordInit["mercure"], a["mercure"], perihelie["mercure"], e["mercure"])
venus = Planet(mass["venus"], coordInit["venus"], a["venus"], perihelie["venus"], e["venus"])
terre = Planet(mass["terre"], coordInit["terre"], a["terre"], perihelie["terre"], e["terre"])
mars = Planet(mass["mars"], coordInit["mars"], a["mars"], perihelie["mars"], e["mars"])
jupiter = Planet(mass["jupiter"], coordInit["jupiter"], a["jupiter"], perihelie["jupiter"], e["jupiter"])
uranus = Planet(mass["uranus"], coordInit["uranus"], a["uranus"], perihelie["uranus"], e["uranus"])
neptune = Planet(mass["neptune"], coordInit["neptune"], a["neptune"], perihelie["neptune"], e["neptune"])