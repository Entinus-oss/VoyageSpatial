import numpy as np
import matplotlib.pyplot as plt
import planetes as pl

G = 6.674*1e-11

class Planet():

    def __init__(self, mass, coord, radius=0):
        self.mass = mass
        self.radius = radius

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

    def distance(self, planet):
        d = np.linalg.norm(self.coord - planet.coord)
        return d
    
    def calcEp(self, planet):
        return -G * self.mass * planet.mass / (self.distance(planet))

    def ep(self, planetArray):
        pot = 0
        for i in range(planetArray.size):
            pot += self.calcEp(planetArray[i])
        return pot

    def ec(self):
        return 1/2 * self.mass * np.linalg.norm(self.speed)**2

    def em(self, planetArray):
        return self.ec() + self.ep(planetArray)  

def a(sondeCoord, planetArray):
    acceleration = 0
    for i in range(planetArray.size):
        acceleration += G * planetArray[i].mass * (planetArray[i].coord - sondeCoord)/ np.linalg.norm(planetArray[i].coord - sondeCoord)**3
    return acceleration

def derivee_u (u, t, planet) :
    # Initialisation de la dérivée
    du = np.empty(u.shape)
    norm = ((planet.x-u[0])**2+(planet.y-u[1])**2)**(1/2)
    # D ́eriv ́ee
    du[0] = u[2]
    du[1] = u[3]
    du[2] = G * planet.mass * (planet.x- u[0])/ norm**3
    du[3] = G * planet.mass * (planet.x - u[1])/ norm**3

    return du

def RK4(u_ini, sonde, planetArray, T, h):

    N = int(T/h)
    t = np.linspace(0, T, N)

    # Initialisation du tableau
    u = np.empty((4, N)) # Condition initiale
    u[:, 0] = u_ini

    for i in range(N-1):

        d1 = derivee_u(u[:, i], t[i], planet)
        d2 = derivee_u(u[:, i] + d1 * h/2, t[i] + h/2, planet)
        d3 = derivee_u(u[:, i] + d2 * h/2, t[i] + h/2, planet)
        d4 = derivee_u(u[:, i] + d3 * h, t[i] + h, planet)
        u[:, i + 1] = u[:, i] + h / 6 * (d1 + 2 * d2 + 2 * d3 + d4)

        for planet in planetArray:
            dist = np.linalg.norm(u[:2, i] - planet.coord)
            if dist <= planet.radius:
                t = np.linspace(0, T, i+1)
                print("Too close from ", planet.name, " : ", dist, "/", planet.radius)
                break

        print(i, "/", N-2, end="\r")
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

        for planet in planetArray:
            dist = np.linalg.norm(u[:2, i] - planet.coord)
            if dist <= planet.radius:
                t = np.linspace(0, T, i+1)
                print("Too close from ", planet.name, " : ", dist, "/", planet.radius)
                break

        print(i, "/", N-2, end="\r")
    return t, u

massPlanet = {"Terre" : 1e24, "Lune" : 7.6e24} #kg
dTerreLune = 3.84400 * 1e4 # m
vMax = 1.7 * 1e4 # km/s


day = 86400 # s
hour = 3600 # s
minute = 60 # s
week = 7 * day 
year = 365 * day

def main():
    theta = 3 * np.pi / 4 # 0 < theta < 2 * np.pi
    vDir = np.array([np.cos(theta), -np.sin(theta)])

    terre = Planet(massPlanet["Terre"], [0, 0])#, pl.terre.radius)
    lune = Planet(massPlanet["Lune"], [dTerreLune, 0])#, 1.7 * 1e6)

    planetArray = np.array([terre])

    R = 35000 * 1e3 #m
    vyIni = np.sqrt(G * terre.mass / R)
    # sonde = Sonde([dTerreLune/2, 0], vMax * vDir, [0, 0])
    sonde = Sonde([R, 0], [0, vyIni], [0, 0])
    # T =  day  # s
    #h = 60 # s
    init = [sonde.x, sonde.y, sonde.vx, sonde.vy]

    #print("x0 :", sonde.coord, "v0 :", sonde.speed)
    #Leapfrog
    #t, u = leapfrog(init, sonde, planetArray, T, h)

    # for i in range(planetArray.size):
    #   plt.plot(planetArray[i].x, planetArray[i].y, 'ro', label=str(i))

    # plt.plot(sonde.x, sonde.y, 'bo', label='sonde start')
    # plt.plot(u[0, -1], u[1, -1], 'b*', label='sonde end')
    # plt.plot(u[0, :], u[1, :])
    '''Calcul de la variation d'énergie maximale sur un plot pour un T 
    (temps d'intégration) donné en fonction du pas h'''

    # N = int(T/h)
    # Ep = np.empty(N)
    # Ec = np.empty(N)
    T = day
    hArray = np.array([30 * minute, hour, day])
    errArray = np.ones(len(hArray))
    # Construction du tableau des temps 
    timeArray = np.arange(1,500) * week 

    print(hArray)
    
    for i in range(len(hArray)) :
        print("\n h =", hArray[i])
        for j in range(len(timeArray)) :
            print("\n t =", timeArray[j])
            t, u = leapfrog(init, sonde, planetArray, timeArray[j], hArray[i])
            N = t.size

            Ep = np.empty(N)
            Ec = np.empty(N)
            Em = np.empty(N)

            for n in range(u.shape[1]) :
                sonde.coord = u[:2, n]
                sonde.speed = u[2:4, n]
                Ep[n] = sonde.ep(planetArray)
                Ec[n] = sonde.ec()
                Em[n] = Ep[n] + Ec[n]
        
            #print(Em)  
            EmMax = np.max(Em) 
            EmMin = np.min(Em)
            errArray[i] = np.abs(EmMax - EmMin)

    # plt.plot(hArray, errArray, 'o', label = r" $ \vert E_{max} - E_{min} \vert $")
        plt.plot(t, Em/Em[0], 'o-', label = "h = " + str(hArray[i]) + " T = " + str(timeArray[j]))
    NArray = np.arange(N)
    # print(len(NArray))
    # print(len(Em))
    #plt.plot(NArray, Ep, label = "énergie potentielle")
    #plt.plot(NArray, Ec, label = "énergie cinétique")
    #plt.plot(NArray, Em, label = 'énergie mécanique') 
    #plt.title("Écart max entre valeurs de l'énergie mécanique pour 200 jours")
    #plt.xlabel("h")
    #plt.ylabel('Ecart maximal entre toutes les valeurs de Em')
    plt.grid()
    plt.legend()
    plt.show()

    # for i in range(u.shape[1]):
    #     sonde.coord = u[:2, i]
    #     sonde.speed = u[2:4, i]
    #     Ep[i] = sonde.ep(planetArray)
    #     Ec[i] = sonde.ec()
    
    # normalizedEp = Ep / np.max(Ep)
    # normalizedEc = Ec / np.max(Ec)

    # plt.plot(t, Ep, label='Ep')
    # plt.plot(t, Ec, label='Ec')
    # Em = Ec + Ep 
    # plt.plot(t, Em, label="Em")
    #plt.ylim(-200000, 200000)
    # plt.xlabel("x [m]")
    # plt.ylabel("y [m]")
    return

if __name__ == "__main__":
    main()