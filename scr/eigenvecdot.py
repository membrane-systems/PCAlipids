import numpy as np
import matplotlib.pyplot as plt
import sys


def main(file1, file2):
	file = open(file1, 'r')

	eigenvecs1 = []

	for line in file:
		line_data = line.split()
		for i in range(len(line_data)):
			line_data[i] = float(line_data[i])
		eigenvecs1.append(line_data)

	file.close()

	file = open(file2, 'r')

	eigenvecs2 = []

	for line in file:
		line_data = line.split()
		for i in range(len(line_data)):
			line_data[i] = float(line_data[i])
		eigenvecs2.append(line_data)

	file.close()

	eigenvecs1 = np.array(eigenvecs1)
	eigenvecs2 = np.array(eigenvecs2)

	matrix = [[0.0] * len(eigenvecs1) for i in range(len(eigenvecs1))]

	for i in range(len(eigenvecs1)):
		for j in range(len(eigenvecs1)):
			matrix[i][j] = abs(np.dot(eigenvecs1[i], eigenvecs2[j]))

	plt.imshow(matrix, interpolation='nearest', origin='lower', cmap='gray_r')
	plt.show()


if __name__ == '__main__':
	args = sys.argv[1:]
	if '-evec' in args:
		main(args[args.index('-evec') + 1], args[args.index('-evec') + 2])
	elif '-h' not in args and '-help' not in args:
		print('Missing parameters, try -h for flags\n')
	else:
		print('-evec <two files with eigenvectors (example: -evec eigenvec1.xvg eigenvec2.xvg)>')
