import numpy
import matplotlib.pyplot as plt
from scipy.stats.kde import gaussian_kde
import numpy as np
import matplotlib.patches as mpatches
import math


def get_data_from_file(file_name):
	file = open(file_name, 'r')
	lines = file.readlines()
	file.close()
	flag = False
	lines = lines[:len(lines) - 1]
	while '' in lines:
		lines.remove()
	data = []
	for line in lines:
		if line.find('@') == -1 and line.find('&') == -1:


			data.append(float(line.split()[0]))
		if line.find('&') != -1:
			break

	return np.array(data)


def PDFs(y, data, label, col, style=''):
	KDEpdf = gaussian_kde(data)
	max_ = np.amax(data)
	min_ = np.amin(data)
	K = KDEpdf(y)
	if style == '':
		p= plt.plot(y, KDEpdf(y), 'r', label = label, color = col) 
	else:
		p= plt.plot(y, KDEpdf(y), 'r', label = label, color = col,linestyle = style)
	# plt.show()
	return max(KDEpdf(y)), p

def main(files,outF):
	if len(files) < 2:
		print("Projection for at least 2 simulations have to be provided.\n\
Run pcalipids.py projdistm -h for help")
		return 0
	handles = []
	data = [get_data_from_file(file) for file in files]
	min_ = 10**10
	max_ = -1 * 10**10
	for seq in data:
		min_ = min(min_, min(seq))
		max_ = max(max_, max(seq))
	y = np.linspace(min_, max_, 101)
	max_t = -1 * 10**10
	colors = ['green','blue','red','cyan','olive','orange','pink','yellow']

	for i in range(len(data)):
		max_y, line1 = PDFs(y, data[i],'%s trajectory' % str(i+1), colors[i%len(colors)])
		max_t = max(max_y, max_t)
		handles.append(line1[0])
	
	plt.title("Distributions in general basis")
	plt.xlabel('PC projection value (A)')
	plt.ylabel('Probability density (a.u.)')
	plt.ylim([0,max_t + 0.05])
	plt.xlim([min_, max_])
	plt.legend(handles = handles, ncol = 2)
	plt.savefig(outF)
	print('Picture saved as '+outF)
