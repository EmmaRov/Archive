import numpy as np 
import matplotlib.pyplot as plt 

part1 = np.loadtxt('part_1_data.txt')
part2 = np.loadtxt('part_2_data.txt')
part3 = np.loadtxt('part_3_data.txt')

# Maybe prettyfy the plot a bit
ax = plt.figure().add_subplot(projection='3d')
ax.plot(part1[:,0], part1[:,1], part1[:,2])
ax.plot(part2[:,0], part2[:,1], part2[:,2])
ax.plot(part3[:,0], part3[:,1], part3[:,2])
plt.show()