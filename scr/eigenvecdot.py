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

	matrix = [[0.0] * 100 for i in range(100)]

	for i in range(100):
		for j in range(100):
			matrix[i][j] = abs(np.dot(eigenvecs1[i], eigenvecs2[j]))


	fig = plt.figure()
	ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
	ax2 = fig.add_axes([0.22, 0.64, 0.24, 0.24])
	ax1.imshow(matrix, interpolation='nearest', origin='lower', cmap='gray_r')
	plt.sca(ax1)
	ax1.set_xticks([0,9,19,29,39,49,59,69,79,89,99]);
	ax1.set_yticks([0,9,19,29,39,49,59,69,79,89,99]);
	ax1.set_xticklabels([1,10,20,30,40,50,60,70,80,90,100]);
	ax1.set_yticklabels([1,10,20,30,40,50,60,70,80,90,100]);
	with open('eigenvecdot.dat', 'w') as file:
		flag = 0
		for i in range(len(matrix)):
			for j in range(len(matrix)):
				file.write(str(matrix[i][j]) + ' ')
				flag += 1
				if flag == 3:
					flag = 0
					file.write('\n')
	print('Wrote eigevectors dot product matrix in "eigenvecdot.dat"')

	matrix = [[0.0] * 10 for i in range(10)]

	for i in range(10):
		for j in range(10):
			matrix[i][j] = abs(np.dot(eigenvecs1[i], eigenvecs2[j]))
	ax2.imshow(matrix, interpolation='nearest', origin='lower', cmap='gray_r')
	plt.sca(ax2)
	ax2.set_xticks([0,4,9]);
	ax2.set_yticks([0,4,9]);
	ax2.set_xticklabels([1,5,10]);
	ax2.set_yticklabels([1,5,10]);
	plt.sca(ax1)
	plt.ylabel('PCALipids eigenvectors for first trajectory')
	plt.xlabel('PCALipids eigenvectors for second trajectory')
	plt.savefig('eigenvecdot.png')
	print('Picture saved in "eigenvecdot.png"')

# if __name__ == '__main__':
# 	args = sys.argv[1:]
# 	if '-evec' in args:
# 		main(args[args.index('-evec') + 1], args[args.index('-evec') + 2])
# 	elif '-h' not in args and '-help' not in args:
# 		print('Missing parameters, try -h for flags\n')
# 	else:
# 		print('-evec <two files with eigenvectors (example: -evec eigenvec1.xvg eigenvec2.xvg)>')
