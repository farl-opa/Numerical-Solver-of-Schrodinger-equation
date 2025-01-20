import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import rc
import os

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Construct the full path to 'mides.txt'
file_path_mides = os.path.join(script_dir, 'mides.txt')

# Read 'mides.txt'
with open(file_path_mides, 'r') as f:
    mides = f.read()
    mides = mides.split('\t')

# Parse T, X, Y from mides
T = int(mides[0])
X = int(mides[1])
Y = int(mides[2])

# Calculate mid-point and other constants
punt_mig = int(X / 2) if X % 2 == 0 else int((X - 1) / 2)
dx = 0.01
L = (X - 1) * dx

# Initialize matrix
matrix = np.zeros((T, X, Y))

# Construct the full path to 'gauss.txt' (adjust if needed)
file_path_gauss = os.path.join(script_dir, 'gaussT.txt')

# Read and parse 'gauss.txt'
with open(file_path_gauss, 'r') as f:
    gauss = f.read()
    gauss = gauss.split('\n\n\n')

for i in range(len(gauss) - 1):
    gauss[i] = gauss[i].split("\n")
    gauss[i].pop(0)  # Remove the first element
    for j in range(len(gauss[i])):
        gauss[i][j] = gauss[i][j].split('\t')
        gauss[i][j].pop()  # Remove the last element
        matrix[i][j] = gauss[i][j]

# Find the maximum value in the matrix
punt_max = np.max(matrix)

# Set up meshgrid for plotting
x = np.linspace(-L / 2, L / 2, X)
y = np.linspace(-L / 2, L / 2, Y)
x, y = np.meshgrid(x, y)

# Create the plot
fig = plt.figure(dpi=100, figsize=(14, 12))
ax = plt.axes(projection='3d')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("$|psi|^2$")
ax.set_zlim3d(0, punt_max)
ax.view_init(90, 90)
ax.ticklabel_format(style="sci")


# Animation function
def animate(i):
    ax.collections.clear()
    ax.plot_surface(x, y, matrix[i], cmap='jet', rstride=1, cstride=1, linewidth=0)
    return ()


# Set up animation
anim = animation.FuncAnimation(fig, animate, frames=T, interval=30, blit=True)
rc('animation', html='jshtml')
plt.close()

file_path_anim = os.path.join(script_dir, 'animation.gif')

# Save animation as GIF
anim.save(file_path_anim, writer='pillow', fps=10)
