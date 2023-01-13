import numpy as np
import matplotlib.pyplot as plt
import planetes as pl

G = 6.674*1e-11

class Planet():

    def __init__(self, name, mass, coord, radius=0):
        self.mass = mass
        self.radius = radius
        self.name = name

        self.coord = coord
        self.x = self.coord[0]
        self.y = self.coord[1]

    
class Sonde():

    def __init__(self, coord, speed, acc, mass=100):
        self.coord = np.array(coord)
        self.speed = np.array(speed)
        self.acc = np.array(acc)

        self.x = coord[0]
        self.y = coord[1]       
        self.vx = speed[0]
        self.vy = speed[1]

        self.mass = mass
        self.potentials = np.array([])
        self.activePotentials = np.array([])
        self.emTrack = []

    def distance(self, planet):
        return np.linalg.norm(self.coord - planet.coord)
    
    def calcEp(self, planet):
        """Calcule Ep"""
        return -G * self.mass * planet.mass / (self.distance(planet))

    def ep(self, planetArray):
        """Calcule Ep, somme sur l'ensemble des éléments de planetArray"""
        pot = 0
        for i in range(planetArray.size):
            pot += self.calcEp(planetArray[i])
        return pot

    def ec(self):
        """Calcule Ec"""
        return 1/2 * self.mass * np.linalg.norm(self.speed)**2

    def em(self, planetArray):
        """Calcule Em"""
        return self.ec() + self.ep(planetArray)  


def a(sondeCoord, planetArray):

    acceleration = 0
    for i in range(planetArray.size):
        acceleration += G * planetArray[i].mass * (planetArray[i].coord - sondeCoord)/ np.linalg.norm(planetArray[i].coord - sondeCoord)**3
    return acceleration

def leapfrog(u_ini, sonde, planetArray, T, h):
    N = int(T/h)
    t = np.linspace(0, T, N)

    # Initialisation du tableau
    u = np.empty((4, N)) # Condition initiale
    u[:, 0] = u_ini    

    # Coefficient d'intégration
    w0 = - 2**(1/3) / (2 - 2**(1/3))
    w1 = 1 / (2 - 2**(1/3))

    c = np.array([w1 / 2, (w0 + w1) / 2, (w0 + w1) / 2, w1/2])
    d = np.array([w1, w0, w1])

    e = sonde.em(planetArray)
    sonde.emTrack.append(e)

    for i in range(N-1):

        #4th order Yoshida integrator
        
        r1 = u[:2, i] + c[0] * u[2:4, i] * h
        v1 = u[2:4, i] + d[0] * a(r1, planetArray) * h

        r2 = r1 + c[1] * v1 * h
        v2 = v1 + d[1] * a(r2, planetArray) * h

        r3 = r2 + c[2] * v2 * h
        v3 = v2 + d[2] * a(r3, planetArray) * h

        u[:2, i+1] = r3 + c[3] * v3 * h
        u[2:4, i+1] = v3 

        sonde.coord = u[:2, i+1]
        sonde.speed = u[2:4, i+1]

        e = sonde.em(planetArray)
        sonde.emTrack.append(e)

        print(i, "/", N-2, end="\r")

        if tooClose(i, u, planetArray):
            #t = np.linspace(0, T, i+1)
            print("Too close !")
            #break

        # if np.abs(np.max(sonde.emTrack) - np.min(sonde.emTrack)) >= 1:
        #     print("Warning : divergence of EM value")
        #     return t[:i], u[:, :i]




    return t, u

def tooClose(i, u, planetArray):
    for planet in planetArray:
        dist = np.linalg.norm(u[:2, i] - planet.coord)
        if dist <= planet.radius:
            return True
            

massPlanet = {"Terre" : 1e24, "Lune" : 7.6e24} #kg
dTerreLune = 3.84400 * 1e4 # m
vMax = 1.7 * 1e4 # km/s


day = 86400 # s
hour = 3600 # s
minute = 60 # s
year = 365 * day

def main():

    terre = Planet("terre", massPlanet["Terre"], [0, 0], pl.terre.radius)
    #lune = Planet(massPlanet["Lune"], [dTerreLune, 0])#, 1.7 * 1e6)

    planetArray = np.array([terre])

    R = 35000 * 1e3 #m
    vyIni = np.sqrt(G * terre.mass / R)
    # sonde = Sonde([dTerreLune/2, 0], vMax * vDir, [0, 0])

    theta = np.pi / 2 # 0 < theta < 2 * np.pi
    vDir = np.array([np.cos(theta), -np.sin(theta)])

    sonde = Sonde([R, 0], [0, vyIni], [0, 0])
    T = 50 * year# s
    h = hour # s
    N = int(T/h)
    init = [sonde.x, sonde.y, sonde.vx, sonde.vy]

    print("x0 :", sonde.coord, "v0 :", sonde.speed)

    #RK4
    #t, u = RK4(init, derivee_u, terre, T, h)

    #Leapfrog
    # t, u = leapfrog(init, sonde, planetArray, T, h)

    # #Plotting...
    # fig, ax = plt.subplots(figsize=(5, 5))

    # r = np.linspace(0, 2 * np.pi, 1000)

    # for i in range(planetArray.size):
    #     ax.plot(planetArray[i].x, planetArray[i].y, 'ro', label=str(i))
    #     ax.plot(pl.terre.radius * np.cos(r), pl.terre.radius * np.sin(r), '--')

    # ax.plot(sonde.x, sonde.y, 'bo', label='sonde start')
    # ax.plot(u[0, -1], u[1, -1], 'b*', label='sonde end')
    # ax.plot(u[0, :], u[1, :])
    # ax.axis('equal')
    # ax.set_xlabel("x [m]")
    # ax.set_ylabel("y [m]")
    # plt.show()

    # Em = np.empty(t.size)

    # for i in range(u.shape[1]):
    #     sonde.coord = u[:2, i]
    #     sonde.speed = u[2:4, i]
    #     Em[i] = sonde.em(planetArray)

    # diffRelatEm = 2 * np.abs(np.max(Em) - np.min(Em)) / np.abs(np.max(Em) + np.min(Em))
    # print("Em_max / Em_min", diffRelatEm)
    # plt.plot(t, Em, label="Em")

    # plt.legend()
    # plt.grid()
    # plt.show()

    return 

if __name__ == "__main__":
    main()