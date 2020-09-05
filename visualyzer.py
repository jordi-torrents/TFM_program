import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

directory='results/VicsekM/configuration/eta15/'
files = os.listdir(directory)
L=np.max(pd.read_csv(directory+files[0], sep='\s+', header=None).T.values)

fig,ax = plt.subplots(1,figsize=(5,5), dpi=150)
ax.set(xlim=(0, L), ylim=(0, L), aspect=1, xticks=(), yticks=())

line, = ax.plot([], [], 'k.', ms=1)
plt.tight_layout()
def animate(i):
    line.set_data(pd.read_csv(directory+files[i], sep='\s+', header=None).T.values)
    return line,

ani = animation.FuncAnimation(fig, animate, frames=5000, interval=13, blit=True, repeat=True)
plt.show()
# ani.save('all_modes_rho40.mp4', fps=25)
