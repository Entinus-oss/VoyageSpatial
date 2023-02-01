import numpy as np
def leapfrog(u_ini, planetArray, t, h):
    N = t.size
    # Initialisation du tableau
    u = np.empty((4, 1)) # Condition initiale
    u[:, 0] = u_ini    

    # Coefficient d'intégration
    w0 = - 2**(1/3) / (2 - 2**(1/3))
    w1 = 1 / (2 - 2**(1/3))

    c = np.array([w1 / 2, (w0 + w1) / 2, (w0 + w1) / 2, w1/2])
    d = np.array([w1, w0, w1])

    em = np.zeros(t.size)
    for planet in planetArray:
        em[0] += calculateEm(u[:2, 0], u[2:4, 0], planet.orbit[:, 0])

    tTracker = 0
    i = 0
    T = t[-1]

    while tTracker < T:

        h = variableTimestep(i, h, u[:2, i], planetArray, week, hour)

        #4th order Yoshida integrator
        
        r1 = u[:2, i] + c[0] * u[2:4, i] * h

        a1 = 0
        for planet in planetArray:
            a1 += a(r1, planet.orbit[:, i])

        v1 = u[2:4, i] - d[0] * a1 * h

        r2 = r1 + c[1] * v1 * h

        a2 = 0
        for planet in planetArray:
            a2 += a(r2, planet.orbit[:, i])

        v2 = v1 - d[1] * a2 * h

        r3 = r2 + c[2] * v2 * h

        a3 = 0
        for planet in planetArray:
            a3 += a(r3, planet.orbit[:, i])

        v3 = v2 - d[2] * a3 * h
        r3 = r3.flatten()
        v3 = v2.flatten()

        x = r3[0] + c[3] * v3[0] * h
        y = r3[1] + c[3] * v3[1] * h
        vx = v3[0]
        vy = v3[1]
        
        u_new = np.array([[x], [y], [vx], [vy]])

        #print("\nu", u[:,], "\nu_new", u_new, "\nappend", np.append(u, u_new, axis=1))
        u = np.append(u, u_new, axis=1)

        for planet in planetArray:
            em[i+1] += calculateEm(u[:2, i+1], u[2:4, i+1], planet.orbit[:, i+1])
    
        # print(i, "/", N-2, end="\r")
        print(tTracker, "/", T, end="\r")
        if np.linalg.norm(planetCoord - u[:2, i+1]) <= MAX_RADIUS:
            print("Error : collision with planet's surface")
            u[:2, i+1] = u[:2, i]
            return t, u, em
        
        i += 1
        tTracker += h

    print("\n")
    return t, u, em

a = np.array([[1,2,3]])
a = np.append(a[0], 3)

print(a)