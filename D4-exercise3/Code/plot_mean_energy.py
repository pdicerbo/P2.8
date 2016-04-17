import matplotlib.pyplot as plt
import numpy as np
import os

all_files    = os.listdir(".")
# all_files    = os.listdir() # uncomment this line if you use python 3.x

energy_files = []

for namef in all_files:
    if namef[0:8] == 'energies':
        energy_files.append(namef)

N = len(energy_files)

t_min = 1.
t_max = 1.5
temperature = np.zeros(N)
mean_energy = np.zeros(N)

for k in range(N):
    temperature[k] = t_min * np.exp(k * np.log(t_max / t_min) / (N-1))

for namef in energy_files:
    # print("Working on", namef)
    print "Working on", namef
    index = int(namef[8:-4])
    data = np.loadtxt(namef)
    energy = data[:,3]
    MyEnergy = np.add.reduce(energy) / len(energy)
    mean_energy[index] = MyEnergy
    # print("mean energy: ", MyEnergy)

# print(mean_energy)
        
plt.figure()
plt.plot(temperature, mean_energy, "-o")
plt.xlabel("T")
plt.ylabel("Mean Energy")
plt.xlim((.8,1.6))

# plt.show()
plt.savefig("mean_energy.png")
