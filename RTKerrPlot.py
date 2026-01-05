import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

print("Loading data...")

with open("shadow_fixed.dat") as f:
    line = f.readline().split()
    W, H = int(line[0]), int(line[1])

data = pd.read_csv("shadow_fixed.dat", sep=" ", skiprows=1, names=["x", "y", "mask", "pattern"])
image = np.zeros((H, W))

for index, row in data.iterrows():
    if row['mask'] == 1:
        image[int(row['y']), int(row['x'])] = row['pattern']
    else:
        image[int(row['y']), int(row['x'])] = 0.0 # Black Hole

plt.figure(figsize=(10, 6), facecolor='black')
plt.imshow(image, cmap='inferno', origin='lower', extent=[-9, 9, -4.5, 4.5])
plt.title(f"Ray Traced Kerr Black Hole (Spin a=0.999)\nGeodesic Integration (RK4)", color='white')
plt.xlabel("Impact Parameter $\\alpha$", color='white')
plt.ylabel("Impact Parameter $\\beta$", color='white')
plt.tick_params(colors='white')
plt.grid(False)
plt.savefig("BH.png", dpi=300, bbox_inches='tight')
plt.show()