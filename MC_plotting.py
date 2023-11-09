import numpy as np 
import matplotlib.pyplot as plt 

Energy_data = np.loadtxt("Energy_samples.txt")
Mag_data = np.loadtxt("Mag_samples.txt")
print(np.shape(Energy_data))

N = len(Energy_data)

n = np.arange(N)

plt.plot(n, Energy_data)
plt.plot(n, Mag_data)
plt.show()

plt.hist(Energy_data, N)
plt.show()

plt.hist(Mag_data, N)
plt.show()