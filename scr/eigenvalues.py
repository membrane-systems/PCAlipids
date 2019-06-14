import numpy as np
import matplotlib.pyplot as plt
import sys

# Load eigenvalues
def get_data(file_name):
	file = open(file_name, 'r')
	lines = file.readlines()
	eigenvals = []
	for i in range(1, len(lines) - 1):
		eigenvals.append(float(lines[i].split()[1]))
	return eigenvals

# Plot eigenvalues
def eig_pic(eigenvals,outF):
	plt.loglog(np.linspace(1, len(eigenvals), len(eigenvals)), eigenvals)
	plt.xlim([1,len(eigenvals)])
	plt.ylabel(r'Eigenvalues ($\AA^2$)')
	plt.xlabel("Component #")
	plt.ylim([10 ** -5, 1])
	plt.savefig(outF)

# Plot cumulative eigenvalues
def eig_pic_cumulative(eigenvals,outF):
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
	plt.ylabel("Cumulative")
	plt.xlabel("Component #")
	plt.yticks([0.1, 0.3, 0.5, 0.7, 0.9, 1.0])
	plt.ylim([0.0, 1.03])
	plt.xlim([1, len(eigcum)])
	plt.savefig(outF)


def main(file_name, cumulative,outF):
	if not file_name:
		print("Eigenvalues have to be provided.\n\
Run pcalipids.py evals -h for help")
		return 0

	eigenvals = get_data(file_name)
	if cumulative:
		eig_pic_cumulative(eigenvals,outF)
	else:
		eig_pic(eigenvals,outF)

