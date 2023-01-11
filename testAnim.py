import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = ax.plot([], [], 'ro')

def init():
    ax.set_xlim(0, 2*np.pi)
    ax.set_ylim(-1, 1)
    return ln,

def update(frame):
    print(frame)
    ydata.append(frame)
    ln.set_data(xdata, ydata)
    return ln,
    
T, N = 1e9, 10000
ani = FuncAnimation(fig, update, frames=np.arange(0, T, T/N, dtype=np.int),
                    init_func=init, blit=True, interval=17)
plt.show()