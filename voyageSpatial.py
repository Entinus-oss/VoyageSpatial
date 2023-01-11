import numpy as np
import matplotlib.pyplot as plt

G = 6.674*1e-11

class Planet():

    def __init__(self, m, coord, v=0, radius=0):
        self.m = m
        self.radius = np.array(radius)
        self.coord = np.array(coord)
        self.x = coord[0]
        self.y = coord[1]
    
class Sonde():

    def __init__(self, coord, speed, acc, mass=100):
        self.coord = np.array(coord)
        self.speed = np.array(speed)
        self.acc = np.array(acc)

        self.x = coord[0]
        self.y = coord[1]       
        self.vx = speed[0]
        self.vy = speed[1]

        self.m = mass
        self.potentials = np.array([])
        self.activePotentials = np.array([])

    def distance(self, planet):
        d = np.linalg.norm(self.coord - planet.coord)
        return d

    def ep(self, planet):
        pot = G * self.m * planet.m / (self.distance(planet))
        np.append(self.potentials, pot)
        return pot

    def ec(self):
        return 1/2 * self.m * np.linalg.norm(self.speed)**2

    def filt_p(self, potentials):
        return

def a(sondeCoord, planetArray):
    acceleration = 0
    for i in range(planetArray.size):
        acceleration += G * planetArray[i].m * (planetArray[i].coord - sondeCoord)/ np.linalg.norm(planetArray[i].coord - sondeCoord)**3
    return acceleration

def derivee_u (u, t, planet) :
    # Initialisation de la dérivée
    du = np.empty(u.shape)
    norm = ((planet.x-u[0])**2+(planet.y-u[1])**2)**(1/2)
    # D ́eriv ́ee
    du[0] = u[2]
    du[1] = u[3]
    du[2] = G * planet.m * (planet.x- u[0])/ norm**3
    du[3] = G * planet.m * (planet.x - u[1])/ norm**3

    return du

def RK4(u_ini, derivee, planet, T, h):

    N = int(T/h)
    t = np.linspace(0, T, N)

    # Initialisation du tableau
    u = np.empty((4, N)) # Condition initiale
    u[:, 0] = u_ini

    for i in range(N - 1):

        d1 = derivee(u[:, i], t[i], planet)
        d2 = derivee(u[:, i] + d1 * h/2, t[i] + h/2, planet)
        d3 = derivee(u[:, i] + d2 * h/2, t[i] + h/2, planet)
        d4 = derivee(u[:, i] + d3 * h, t[i] + h, planet)
        u[:, i + 1] = u[:, i] + h / 6 * (d1 + 2 * d2 + 2 * d3 + d4)

        if np.linalg.norm(u[:2, i] - planet.coord) < planet.radius:
            t = np.linspace(0, T, i+1)
            break

    return t, u


def leapfrog(u_ini, sonde, planetArray, T, h):
    N = int(T/h)
    t = np.linspace(0, T, N)

    # Initialisation du tableau
    u = np.empty((4, N)) # Condition initiale
    u[:, 0] = u_ini    

    w0 = - 2**(1/3) / (2 - 2**(1/3))
    w1 = 1 / (2 - 2**(1/3))

    c = np.array([w1 / 2, (w0 + w1) / 2, (w0 + w1) / 2, w1/2])
    d = np.array([w1, w0, w1])

    #a = lambda r: G * planet.m * (planet.coord - r)/ np.linalg.norm(planet.coord - r)**3

    for i in range(N - 1):

        #4th order Yoshida integrator
        
        r1 = u[:2, i] + c[0] * u[2:4, i] * h
        v1 = u[2:4, i] + d[0] * a(r1, planetArray) * h

        r2 = r1 + c[1] * v1 * h
        v2 = v1 + d[1] * a(r2, planetArray) * h

        r3 = r2 + c[2] * v2 * h
        v3 = v2 + d[2] * a(r3, planetArray) * h

        u[:2, i+1] = r3 + c[3] * v3 * h
        u[2:4, i+1] = v3 

    return t, u

massPlanet = {"Terre" : 1e24, "Lune" : 7.6e24} #kg
dTerreLune = 384400 #km

def main():

    terre = Planet(massPlanet["Terre"], [0, 0])
    lune = Planet(massPlanet["Lune"], [dTerreLune, 0])

    planetArray = np.array([terre, lune])

    R = 35000 #km
    vyIni = np.sqrt(G * terre.m / R)
    sonde = Sonde([dTerreLune/2, 0], [0, vyIni], [0, 0])
    # sonde = Sonde([100, 100], [-1, 1], [0, 0])
    T = 30
    h = 1e-3
    init = [sonde.x, sonde.y, sonde.vx, sonde.vy]

    #RK4
    #t, u = RK4(init, derivee_u, terre, T, h)

    #Leapfrog
    t, u = leapfrog(init, sonde, planetArray, T, h)

    for i in range(planetArray.size):
      plt.plot(planetArray[i].x, planetArray[i].y, 'ro', label=str(i))

    plt.plot(sonde.x, sonde.y, 'bo', label='sonde start')
    plt.plot(u[0, -1], u[1, -1], 'b*', label='sonde end')
    plt.plot(u[0, :], u[1, :])

    # Ep = np.empty(N)
    # Ec = np.empty(N)

    # for i in range(u.shape[1]):
    #     sonde.coord = u[:2, i]
    #     sonde.speed = u[2:4, i]
    #     Ep[i] = sonde.ep(terre)
    #     Ec[i] = sonde.ec()
    
    # normalizedEp = Ep / np.max(Ep)
    # normalizedEc = Ec / np.max(Ec)

    # # plt.plot(t, normalizedEp, label='Ep')
    # # plt.plot(t, normalizedEc, label='Ec')
    # Em = Ec + Ep 
    # plt.plot(t, Em, label="Em")
    plt.ylim(-200000, 200000)

    plt.legend()
    # plt.xlabel("x")
    # plt.ylabel("y")
    plt.show()
    return

if __name__ == "__main__":
    main()