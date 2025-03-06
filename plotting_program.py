import numpy as np
import matplotlib.pyplot as plt
import glob
import os

file_list = sorted(glob.glob("klein_gordon_output_*.dat"), key=lambda x: int(x.split("_")[-1].split(".")[0]))

if not file_list:
    print("no data has been generated")
    exit()

def load_data(filename):
    data = np.loadtxt(filename)
    return data[:, 0], data[:, 1]

plt.figure(figsize=(10, 5))
for file in file_list:
    x, phi = load_data(file)
    step = file.split("_")[-1].split(".")[0]
    plt.plot(x, phi, alpha=0.5, label=f"Step {step}")

plt.xlabel("$x$")
plt.ylabel("φ(x,t)")
plt.title("Klein-Gordon Wave Propagation")
plt.legend(fontsize=8, loc="upper right", ncol=2)
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 5))

x_init, phi_init = load_data(file_list[0])
plt.plot(x_init, phi_init, label="Step 0", linewidth=2)

x_final, phi_final = load_data(file_list[-1])
plt.plot(x_final, phi_final, label=f"Step {file_list[-1].split('_')[-1].split('.')[0]}", linewidth=2)

plt.xlabel("$x$")
plt.ylabel("φ(x, t)")
plt.title("Klein-Gordon Wave Propagation (Initial & Final Steps)")
plt.legend()
plt.grid(True)
plt.show()
