import numpy as np
import matplotlib.pyplot as plt
import glob

file_list = sorted(glob.glob("klein_gordon_output_*.dat"), key=lambda x: int(x.split("_")[-1].split(".")[0]))

if not file_list:
    print("no data has been generated")
    exit()

time_steps = []
phi_values = []

for file in file_list:
    step = int(file.split("_")[-1].split(".")[0])
    data = np.loadtxt(file)
    x_values = data[:, 0]  
    phi_values.append(data[:, 1]) 
    time_steps.append(step)

phi_values = np.array(phi_values)
time_values = np.array(time_steps)

plt.figure(figsize=(10, 6))
plt.imshow(phi_values, aspect='auto', extent=[x_values[0], x_values[-1], time_values[0], time_values[-1]], cmap='viridis', origin='lower')
plt.colorbar(label="φ(x, t)")
plt.xlabel("$x$")
plt.ylabel("$t$")
plt.title("Time Evolution of φ(x, t)")

plt.savefig("klein_gordon_heatmap.png")
plt.show()
