import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib import cm
from tqdm import tqdm

X = 40
Y = 30
T = 10000

dx = 1
dy = 1
dt = 0.0001

rho = 1.21
c = 343
K = rho*c**2

Q = 0.5 + 0.5 * np.cos(np.arange(-np.pi, np.pi, 2*np.pi/200))


P = np.zeros((X, Y), "float64")
Ux = np.zeros((X+1, Y), "float64")
Uy = np.zeros((X, Y+1), "float64")

P_all = np.zeros((T, X, Y), "float64")
mic = []
Z = rho * c

for n in tqdm(range(T)):
    if n < len(Q):
        P[20, 15] += Q[n]
    P_all[n] = P

    # equation of motion
    Ux[1:X, :] = Ux[1:X, :] - dt/(rho*dx) * (P[1:X, :] - P[:X-1, :])
    Uy[:, 1:Y] = Uy[:, 1:Y] - dt/(rho*dy) * (P[:, 1:Y] - P[:, :Y-1])
    
    # impedance
    Ux[0, :] = P[0, :] / -Z
    Ux[-1, :] = P[-1, :]/ Z

    Uy[:, 0] = P[:, 0] / -Z
    Uy[:, -1] = P[:, -1]/ Z

    # equation of continuity
    P[:X, :Y] = P[:X, :Y] - K*dt/dx * (Ux[1:X+1, :] - Ux[:X, :]) - K*dt/dy * (Uy[:, 1:Y+1] - Uy[:, :Y])
    

# prepare figure
fig = plt.figure()
fig.set_dpi(100)
ax = Axes3D(fig)

x = np.arange(0, dx*X, dx)
y = np.arange(0, dy*Y, dy)
xx, yy = np.meshgrid(y, x)

def animate(i):
   ax.clear()
   ax.plot_surface(xx, yy, P_all[i], rstride=1, cstride=1, cmap=plt.cm.coolwarm, vmax=1, vmin=-1)
   ax.set_zlim(-2, 2)

anim = animation.FuncAnimation(fig, animate, frames=len(P_all) - 1, interval=1, repeat=False)
anim.save("new.mp4", fps = 60)
plt.show()
