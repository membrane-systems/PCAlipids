import numpy as np
import matplotlib.pyplot as plt

file = open('eigenvec_aa_to_map.xvg', 'r')

eigenvecs = []

for line in file:
	line_data = line.split()
	for i in range(len(line_data)):
		line_data[i] = float(line_data[i])
	eigenvecs.append(line_data)

file.close()

eigenvecs = np.array(eigenvecs)

matrix = [[0.0] * len(eigenvecs) for i in range(len(eigenvecs))]

for i in range(len(eigenvecs)):
	for j in range(len(eigenvecs)):
		matrix[i][j] = np.dot(eigenvecs[i], eigenvecs[j])

plt.imshow(matrix, interpolation='nearest', origin='lower')
plt.show()