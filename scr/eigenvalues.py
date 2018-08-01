import numpy as np
import matplotlib.pyplot as plt
import sys


def get_data(file_name):
	file = open(file_name, 'r')
	lines = file.readlines()
	eigenvals = []
	for i in range(1, len(lines) - 1):
		print(lines[i].split()[1])
		eigenvals.append(float(lines[i].split()[1]))
	return eigenvals


def eig_pic(eigenvals):
	plt.loglog(np.linspace(1, len(eigenvals), len(eigenvals)), eigenvals)
	plt.xlim([1,len(eigenvals)])
	plt.ylabel("Eigenvalues nm^2")
	plt.xlabel("Component")
	plt.ylim([10 ** -4, 10])
	plt.show()


def eig_pic_cumulative(eigenvals):
	eigsum = sum(eigenvals)
	eigcum = [eigenvals[0]]
	for i in range(1,len(eigenvals)):
		eigcum.append(eigcum[i - 1] + eigenvals[i])
	for i in range(len(eigcum)):
		eigcum[i] /= eigsum
	plt.semilogx(np.linspace(1, len(eigcum), len(eigcum)), eigcum)
	plt.axhline(y = 0.5, color = 'black', linestyle = '--')
	plt.axhline(y = 0.9, color = 'black', linestyle = '--')
	plt.axhline(y = 1, color = 'black', linestyle = '--')
	plt.ylabel("Cumulative nm^2")
	plt.xlabel("Component")
	plt.yticks([0.1, 0.3, 0.5, 0.7, 0.9, 1.0])
	plt.ylim([0.0, 1.03])
	plt.xlim([1, len(eigcum)])
	plt.show()


def main(file_name, cumulative = False):
	eigenvals = get_data(file_name)
	if cumulative:
		eig_pic_cumulative(eigenvals)
	else:
		eig_pic(eigenvals)


# if __name__ == '__main__':
# 	args = sys.argv[1:]
# 	if len(args) == 2:
# 		main(args[0], args[1])
# 	else:
# 		main(args[0])
