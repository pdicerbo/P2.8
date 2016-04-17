import numpy as np
import matplotlib.pyplot as plt

file1 = "../D1-exercise1/Code/timing_def.dat" # first implementation
file2 = "../D2-exercise2/Code/timing_def.dat" # with linked cells

mat1 = np.loadtxt(file1)
mat2 = np.loadtxt(file2)

npes  = mat1[:,0]
time1 = mat1[:,1]
time2 = mat2[:,1]

plt.figure()
plt.plot(npes, time1, "-o", label = "naive implementation")
plt.plot(npes, time2, "-o", label = "with linked cells")
plt.xlabel("Number of processes")
plt.ylabel("Time (s)")
plt.title("Timing (natoms = 62500)")
plt.legend()
plt.savefig("timing.png")

plt.close("all")

speed1 = time1[0] / time1
speed2 = time2[0] / time2

plt.plot(npes, npes, label = "ideal scaling")
plt.plot(npes, speed1, "-o", label = "naive implementation")
plt.plot(npes, speed2, "-o", label = "with linked cells")
plt.xlabel("Number of processes")
plt.ylabel("Speedup")
plt.title("Scaling (natoms = 62500)")
plt.legend()

plt.savefig("speedup.png")
