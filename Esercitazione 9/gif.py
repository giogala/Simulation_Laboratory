
from matplotlib.animation import PillowWriter, FuncAnimation
import os
import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
plt.rcParams.update({
    'font.family':'cmr10',
    'mathtext.fontset': 'cm',
    'axes.formatter.use_mathtext': True,
    'figure.figsize': [6.0, 6.0],
    'axes.labelsize': 20,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'font.size':10,
    'savefig.directory':'./Images'
})
gen = 300
frames = []
output_dir = 'CIRCLE/BEST'
os.makedirs(output_dir, exist_ok=True)
for j in range(gen):
    t=np.loadtxt("CIRCLE/PATH/element_"+str(j)+"_0.txt",delimiter='\t',skiprows=0)
    x = t[:,1]
    y = t[:,2]
    l = t[:,0]
    f = np.arange(len(x))
    fig, ax = plt.subplots()
    plt.scatter(x,y,c=f,cmap='plasma',zorder=2)
    plt.plot(x,y,color='grey',marker='o',linewidth=0.8,zorder=1)
    circle1 = plt.Circle((0,0),1,color='r',alpha=0.5,fill=False,zorder=0)
    for i,txt in enumerate(l):
        plt.annotate(int(txt),(x[i]*1.05-0.02,y[i]*1.05-0.02))
    ax.add_patch(circle1)
    frame_filename = os.path.join(output_dir, f'frame_{j:03d}.png')
    plt.savefig(frame_filename)
    frames.append(frame_filename)
    plt.close()
    print(j)

fig, ax = plt.subplots(figsize=(6, 6))
ax.axis('off')  # Nasconde gli assi anche per l'animazione finale
im = ax.imshow(plt.imread(frames[0]))


# Funzione di aggiornamento per ogni frame
def update(frame):
    img = plt.imread(frames[frame])
    im.set_data(img)
    return im,

# Crea l'animazione
ani = FuncAnimation(fig, update, frames=len(frames), interval=100, blit=True)

# Salva l'animazione come GIF
gif_filename = 'circle.gif'
ani.save(output_dir+"/"+gif_filename, writer=PillowWriter(fps=10))

# Rimuovi i frame temporanei
for frame in frames:
    os.remove(frame)