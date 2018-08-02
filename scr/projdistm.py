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
	return p

def main(file1, file2):
	handles = []
	data1 = get_data_from_file(file1)
	data2 = get_data_from_file(file2)
	min_ = min([min(data1), min(data2)])
	max_ = max([max(data1), max(data2)])
	y = np.linspace(min_, max_, 101)
	line1, = PDFs(y, get_data_from_file(file1),'First trajectory', 'blue')
	handles.append(line1)
	
	
	line1, = PDFs(y, get_data_from_file(file2),'Second trajectory', 'green')
	handles.append(line1)
	#line1, = PDFs(y, get_data_from_file('projection_0004_map_12_or.xvg'), 'original mapping','black', '--')
	#handles.append(line1)
	#PDFs(y, get_data_from_file('projection_0008_map_16.xvg'), 'black', '--')
	#line1, = PDFs(y, get_data_from_file('projection_slipids.xvg'),'Slipids', (52./255,167./255,7./255))
	#handles.append(line1)
	plt.title("Distributions in general basis")
	plt.xlabel('PC projection value (A)')
	plt.ylabel('Probability density (a.u.)')
	plt.ylim([0,0.45])
	plt.xlim([min_, max_])
	plt.legend(handles = handles, ncol = 2)
	plt.savefig('Distributions_in_general_basis.png')
	print('Picture saved as "Distributions_in_general_basis.png"')