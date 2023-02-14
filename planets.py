import numpy as np
import matplotlib.pyplot as plt

Ms = 1.989 * 1e30 # masse du soleil (kg)
G = 6.674*1e-11 # constante gravitaionnel

class Planet():

    def __init__ (self, mass, orbit, hillSphereRadius, radius=0):
        """
        Classe "Planet" : permet d'accéder facilement aux attributs des planètes.
        """
        self.mass=mass
        self.orbit=orbit
        self.raidus=radius
        self.hillSphereRadius = hillSphereRadius

def calculateInitialPhase(coordInit, a):

        """
        coordInit : np.array([x, y]) coordonnées initiales des planètes.
        a : float (m) demi grand-axe des planètes
        Calcule la phase initiale des planètes du système solaire par rapport à la date de lancer de Voyager 2.
        """
        coordInitNorm = np.linalg.norm(coordInit)
        aNorm = np.linalg.norm(np.array([a, 0]))
        c = np.linalg.norm(coordInit - np.array([a, 0]))

        y = coordInit[1]

        if y <= 0 :
                phaseInit = 2 * np.pi - np.arccos((coordInitNorm**2 + aNorm**2 - c**2)/(2 * coordInitNorm * aNorm))        
        else : 
                phaseInit = np.arccos((coordInitNorm**2 + aNorm**2 - c**2)/(2 * coordInitNorm * aNorm)) 

        return phaseInit    

def calculateOrbit(time, phaseInit, a, perihelie, e, mass):

        """
        time : np.array([0, ..., T]) timeline sur laquelle l'orbite des planètes est calculé
        phaseInit : float phase initiale des planètes
        a : float (m) demi grand-axe des planètes
        périhélie : flaot (m) périhélie des planètes
        e : flaot (m) eccentricité des planètes
        mass : flaot (kg) masse des planètes
        """
        u = - a + perihelie
        b = a * np.sqrt(1 - e**2)
        period = 2 * np.pi * np.sqrt(a**3 / (Ms + mass) / G)
        omega = 2 * np.pi / period

        orbit = [u + a * np.cos(omega * time + phaseInit), b * np.sin(omega * time + phaseInit)]

        return np.array(orbit).T

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

hillSphereRadius = {"mercure" : a["mercure"] * (1-e["mercure"]) * (mass["mercure"] / Ms / 3)**(1/3),
                "venus" : a["venus"] * (1-e["venus"]) * (mass["venus"] / Ms / 3)**(1/3),
                "terre" : a["terre"] * (1-e["terre"]) * (mass["terre"] / Ms / 3)**(1/3),
                "mars" : a["mars"] * (1-e["mars"]) * (mass["mars"] / Ms / 3)**(1/3),
                "jupiter" : a["jupiter"] * (1-e["jupiter"]) * (mass["jupiter"] / Ms / 3)**(1/3),
                "saturne" : a["saturne"] * (1-e["saturne"]) * (mass["saturne"] / Ms / 3)**(1/3),
                "uranus" : a["uranus"] * (1-e["uranus"]) * (mass["uranus"] / Ms / 3)**(1/3),
                "neptune" : a["neptune"] * (1-e["neptune"]) * (mass["neptune"] / Ms / 3)**(1/3)}

perihelie = {"mercure" : 46e9, # m
            "venus" : 107.480e9,
            "terre" : 147.095e9,
            "mars" : 206.650e9,
            "jupiter" : 740.595e9,
            "saturne" : 1357.554e9,
            "uranus" : 2732.696e9,
            "neptune" : 4471.050e9}

coordInit = {"mercure" : [1.9379e10, -6.513e10],
            "venus" : [6.588e10,  8.521612947350608e10],
            "terre" : [1.275e11, -8.269e10],
            "mars" : [1.573e11,  1.511e11],
            "jupiter" : [1.243e11,  7.509e11],
            "saturne" : [-1.066e12,  8.643e11],
            "uranus" : [-2.083e12, -1.842e12],
            "neptune" : [-1.132e12, -4.386e12]}           

coordFinal = {"mercure" : [-5.709349323145684e7, -2.154788582251384e7], # 25 aout 1981
            "venus" : [-4.826049890952200e7, -9.585580109732188e7], # km
            "terre" : [ 1.346993039508645e8, -7.041218300374475e7],
            "mars" : [-8.890131717038222e6,  2.361420877941326e8],
            "jupiter" : [-7.777890754287263e8, -2.408650801885601e8],
            "saturne" : [-1.403979755879619e9, -2.905642873950667e8],
            "uranus" : [-1.426015556248566e9, -2.428069102412113e9],
            "neptune" : [-4.591963134603668e8, -4.503837037318032e9]} 

phaseInit = {"mercure" : calculateInitialPhase(coordInit["mercure"], a["mercure"]), # m
            "venus" : calculateInitialPhase(coordInit["venus"], a["venus"]),
            "terre" : calculateInitialPhase(coordInit["terre"], a["terre"]),
            "mars" : calculateInitialPhase(coordInit["mars"], a["mars"]),
            "jupiter" : calculateInitialPhase(coordInit["jupiter"], a["jupiter"]),
            "saturne" : calculateInitialPhase(coordInit["saturne"], a["saturne"]),
            "uranus" : calculateInitialPhase(coordInit["uranus"], a["uranus"]),
            "neptune" : calculateInitialPhase(coordInit["neptune"], a["neptune"])}

### Initialisation des variables de temps ###

minute = 60 # s
hour = 60 * minute
day = 24 * hour
week = 7 * day
month = 30.5 * day
year = 365 * day

T = 1/2*year #2.649e9 # s
dt = 1 # s

time = np.arange(0, T, dt)

### Calcule des orbites ###
orbit = {"mercure" : calculateOrbit(time, phaseInit["mercure"], a["mercure"], perihelie["mercure"], e["mercure"], mass["mercure"]), # m
        "venus" : calculateOrbit(time, phaseInit["venus"], a["venus"], perihelie["venus"], e["venus"], mass["venus"]),
        "terre" : calculateOrbit(time, phaseInit["terre"], a["terre"], perihelie["terre"], e["terre"], mass["terre"]),
        "mars" : calculateOrbit(time, phaseInit["mars"], a["mars"], perihelie["mars"], e["mars"], mass["mars"]),
        "jupiter" : calculateOrbit(time, phaseInit["jupiter"], a["jupiter"], perihelie["jupiter"], e["jupiter"], mass["jupiter"]),
        "saturne" : calculateOrbit(time, phaseInit["saturne"], a["saturne"], perihelie["saturne"], e["saturne"], mass["saturne"]),
        "uranus" : calculateOrbit(time, phaseInit["uranus"], a["uranus"], perihelie["uranus"], e["uranus"], mass["uranus"]),
        "neptune" : calculateOrbit(time, phaseInit["neptune"], a["neptune"], perihelie["neptune"], e["neptune"], mass["neptune"])}

### Instance des objets "Planet" ###
soleil = Planet(Ms, [0, 0], 0)
terre = Planet(mass["terre"], orbit["terre"], hillSphereRadius["terre"])
venus = Planet(mass["venus"], orbit["venus"], hillSphereRadius["venus"])
mars = Planet(mass["mars"], orbit["mars"], hillSphereRadius["mars"])

planetArray = np.array([terre, venus, mars])


### Affichage ####
for planet in planetArray:
        plt.plot(planet.orbit[::100, 0], planet.orbit[::100, 1], 'r-', label = "planet")
        plt.plot(planet.orbit[0, 0], planet.orbit[0, 1], 'bo', label = "planet start")
        plt.plot(planet.hillSphereRadius * np.cos(np.linspace(0, 2 * np.pi, 1000)) + planet.orbit[0, 0], planet.hillSphereRadius * np.sin(np.linspace(0, 2 * np.pi, 1000)) + planet.orbit[0, 1], 'r--')

plt.plot(0, 0, 'ro', label="soleil")

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

plt.grid()
plt.show()