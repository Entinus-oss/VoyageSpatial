import numpy as np
import matplotlib.pyplot as plt

# masse (kg)
# longueur (m)
# temps (s)

M = 5.972e24 # masse du Soleil 
Ms = 1.989 * 1e30 # masse de la Terre
TS = 149.598e9 # distance Terre-Soleil
m = 100 # masse de la sonde
G = 6.674*1e-11 # constante gravitaionnel

# temps
minute = 60 # s
hour = 60 * minute
day = 24 * hour
week = 7 * day
month = 30.5 * day
year = 365 * day

def calculateEm(r, v, rPlanet):
    """
    r : np.array (2) distance de la sonde par rapport à l'origine
    v : np.array (2) vitesse de la sonde
    rPlanet : np.array (2) distance de la planète par rapport à l'origine

    Calcule l'énergie mécanique de la sonde.
    """
    return 1/2 * m * np.linalg.norm(v)**2 - G * M * m / np.linalg.norm(rPlanet - r)

def a(sondeOrbit, planetOrbit):
    """
    sondeOrbit : np.array (2) distance de la sonde par rapport à l'origine
    planetOrbit : np.array (2) distance de la planète par rapport à l'origine

    Calcule l'accélération de la sonde.
    """
    return -G * M * (planetOrbit - sondeOrbit)/ np.linalg.norm(planetOrbit - sondeOrbit)**3

def leapfrog(u_ini, planetArray, t, h):
    """
    u_ini : nparray([x, y, vx, vy]) données initiales de la sonde
    planetArray : npArray([planet1, planet2, ...]) Array contenant des objets Planet
    t : np.array([0, .., T]) échelle de temps
    h : float pas de temps

    Résoud l'équation différentielle pour a. Défini aussi une nouvelle échelle des temps t' non linéaire,
    cela est du à l'algorithme de pas de temps variable.
    """
    N = t.size
    # Initialisation du tableau
    u = np.empty((4, 1)) # Condition initiale
    u[:, 0] = u_ini    

    # Coefficient d'intégration de Yoshida
    w0 = - 2**(1/3) / (2 - 2**(1/3))
    w1 = 1 / (2 - 2**(1/3))

    c = np.array([w1 / 2, (w0 + w1) / 2, (w0 + w1) / 2, w1/2])
    d = np.array([w1, w0, w1])

    # Initialisation de l'énergie mécanique
    em = np.array([0])
    for planet in planetArray:
        em[0] += calculateEm(u[:2, 0], u[2:4, 0], planet.orbit[0, :])

    time = np.array([0])
    i = 0
    T = t[-1]
    timeTracker = 0

    while timeTracker <= T:

        # calcule h avec l'algorithme à pas de temps variable
        h = np.floor(variableTimestep(i, h, u[:2, i], planetArray, month, 10*minute, 1.1, 0.7))
        time = np.append(time, h + time[-1])

        #4th order Yoshida integrator
        
        r1 = u[:2, i] + c[0] * u[2:4, i] * h

        a1 = 0
        for planet in planetArray: # additionne l'accélération de chaque planet de planetArray
            a1 += a(r1, planet.orbit[i, :])

        v1 = u[2:4, i] - d[0] * a1 * h

        r2 = r1 + c[1] * v1 * h

        a2 = 0
        for planet in planetArray: # même chose
            a2 += a(r2, planet.orbit[i, :])

        v2 = v1 - d[1] * a2 * h

        r3 = r2 + c[2] * v2 * h

        a3 = 0
        for planet in planetArray: # même chose
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

        # calcule l'énergie mécanique
        emTemp = 0
        for planet in planetArray:
            # print(planet.orbit[i+1, :], u[:2, i+1])
            emTemp += calculateEm(u[:2, i+1], u[2:4, i+1], planet.orbit[i+1, :])
        em = np.append(em, emTemp)
    
        print(timeTracker, "/", T, end="\r")

        # Verifie que la sonde ne soit pas trop proche de la planète
        for planet in planetArray:
            if np.linalg.norm(planet.orbit[i, :] - u[:2, i+1]) <= MAX_RADIUS:
                print("\nError : collision with planet's surface")
                u[:2, i+1] = u[:2, i]
                return t, u, em

        timeTracker += h # ajoute le pas de temps à la timeline
        i += 1

    print("\n")
    return time, u, em

def variableTimestep(i, h, sondeCoord, planetArray, hMax, hMin, upScale, downScale, hillSphereScaleFactor=3):
    """
    i : int indice en cours
    h : float pas de temps précédent
    sondeCoord : np.array([x , y]) pas de temps précédent 
    planetArray : np.array([planet1, planet2, ...]) arrray contenant des objets Planet
    hMax : float seuil max pour le pas de temps 
    hMin : float seuil min pour le pas de temps 
    upScale : float > 1 facteur de grossissement du pas de temps
    downScale = 0 < float < 1 facteur de rapetissement du pas de temps
    hillSphereScaleFactor : int voir le rapport
    """
    newh = h # copy h
    for planet in planetArray: 
        threshold = hillSphereScaleFactor * planet.hillSphereRadius # vérifie si la sonde est dans le threshold

        if np.linalg.norm(planet.orbit[i, :] - sondeCoord) <= threshold: # si oui réduit le pas de temps
            newh *= downScale
            if newh < hMin:
                newh = hMin
        else: # sinon augmente le pas de temps
            newh *= upScale
            if newh > hMax:
                newh = hMax
    return newh
            
MAX_RADIUS = 6.3e6

class Planet():

    """
    Classe "Planet" : permet d'accéder facilement aux attributs des planètes.
    """

    def __init__ (self, mass, orbit, hillSphereRadius, radius=0):
        self.mass=mass
        self.orbit=orbit
        self.raidus=radius
        self.hillSphereRadius = hillSphereRadius

### Conditions initiales ####
rayonSphereHill = TS * (1-0.0167) * (M / Ms / 3)**(1/3)# m

R = 35000 * 1e3 #m
vyIni = np.sqrt(G * M / R)
vitesseSonde = 1700 # m/s
sondeInitCoord = np.array([0, R])# m
sondeInitSpeed = np.array([0,  vyIni])

vitesseTerre = 30000 # m/s5
u_ini = np.append(sondeInitCoord, sondeInitSpeed)

T = 100*year
dt = 1 #10 * minute

t = np.arange(0, T, dt)

### Orbites des planètes ####
# planet_orbit = np.array([np.cos(2 * np.pi * t / P), np.sin(2 * np.pi * t / P)]).T * R
planet_orbit = np.zeros([t.size, 2])
planet_orbit2 = np.array([np.ones(t.size) * 1e11, np.ones(t.size) * 1.5e11]).T

# planet_orbit = np.array([t * vitesseTerre, np.zeros(t.size)]).T
# print(planet_orbit)

### Instance des objets ####
planet = Planet(M, planet_orbit, rayonSphereHill)
planet2 = Planet(M, planet_orbit2, rayonSphereHill)

planetArray = np.array([planet,])

### Algorithme leapfrog ####
time,u,em = leapfrog(u_ini, planetArray, t, dt)


# Affichage Sonde
plt.plot(u[0, :], u[1, :], '-', label = "Trajectoire de la sonde")
plt.plot(u[0, 0], u[1, 0], 'bo', label = "Position initiale de la sonde")
# plt.plot(u[0, -1], u[1, -1], 'b*', label = "sonde end")

# Affichage Planet
for planet in planetArray:
    plt.plot(planet_orbit[::100, 0], planet_orbit[::100, 1], label = "Trajectoire de la Terre")
    plt.plot(planet.orbit[0, 0], planet.orbit[0, 1], 'ro', label = "Terre")
    plt.plot(MAX_RADIUS * np.cos(np.linspace(0, 2 * np.pi, 1000)) + planet.orbit[0, 0], MAX_RADIUS * np.sin(np.linspace(0, 2 * np.pi, 1000)) + planet.orbit[0, 1], 'r--', label="rayon")
    # plt.plot(planet_orbit[-1, 0], planet_orbit[-1, 1], 'r*', label = "Position finale de la planète")


# plt.xlim(-2e9, 2e9)
# plt.ylim(-2e9, 2e9)
plt.xlabel("x [m]")
plt.ylabel("y [m]")

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.legend()
plt.grid()
plt.title("Ralentissement de la sonde par effet de fronde")
plt.show()

# Affichage Em
emRelatif = np.abs(2 * (np.max(em) - np.min(em)) / (np.max(em) + np.min(em)))
# plt.plot(time, em, '-', label="em")
# print("Em =", emRelatif)

# Affichage vitesse de la Sonde 
# print(u[2:4], u[2, :], u[3, :])
vx = u[2, :]
vy = u[3, :]
v = np.sqrt(np.power(vx, 2) + np.power(vy, 2))
# print(v)

plt.plot(time, v, '-')
plt.xlabel("t [s]")
plt.ylabel("v [m/s]")

plt.title("Évolution temporelle de la vitesse sur un an")
plt.grid()
plt.legend()
plt.show()