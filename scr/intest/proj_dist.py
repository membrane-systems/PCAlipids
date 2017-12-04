import numpy
import matplotlib.pyplot as plt
from scipy.stats.kde import gaussian_kde
import time 
import sys
import numpy as np
import matplotlib.patches as mpatches
import math
from scipy import stats


def get_data_from_file(file_name):
	file = open(file_name, 'r')
	data = []
	for line in file:
		if line.find('@') == -1 and line.find('&') == -1:
			print(line.split())

			data.append(float(line.split()[1]) * 10 / math.sqrt(12))
		if line.find('&') != -1:
			break
	file.close()
	return np.array(data)


def PDFs(y, data, label, col, style=''):
	KDEpdf = gaussian_kde(data)
	max_ = np.amax(data)
	min_ = np.amin(data)
	K = KDEpdf(y)
	print(max_)
	if style == '':
		p= plt.plot(y, K * math.sqrt(12), 'r', label = label, color = col) 
	else:
		p= plt.plot(y, K * math.sqrt(12), 'r', label = label, color = col,linestyle = style)
	return p

def main():
	y = np.linspace(-3 , 6 , 101)
	handles = []
	
	line1, = PDFs(y, get_data_from_file('1500_DOPC_projection.xvg'),'AA DOPC', '#008744')
	handles.append(line1)
	line1, = PDFs(y,get_data_from_file('1510_SOPC_projection.xvg'),'AA SOPC', '#ffa700')
	handles.append(line1)
	line1, = PDFs(y, get_data_from_file('1520_OSPC_projection.xvg'),'AA OSPC', '#0057e7')
	handles.append(line1)
	line1, = PDFs(y, get_data_from_file('1530_DSPC_projection.xvg'),'AA DSPC', '#d62d20')
	handles.append(line1)
	line1, = PDFs(y, get_data_from_file('3100_DOPC_projection.xvg'),'CG DOPC', '#008744', '--')
	handles.append(line1)
	line1, = PDFs(y, get_data_from_file('3110_SOPC_projection.xvg'), 'CG SOPC','#ffa700', '--')
	handles.append(line1)
	line1, = PDFs(y, get_data_from_file('3120_OSPC_projection.xvg'),'CG OSPC', '#0057e7', '--')
	handles.append(line1)
	#line1, = PDFs(y, get_data_from_file('projection_0004_map_12_or.xvg'), 'original mapping','black', '--')
	#handles.append(line1)
	#PDFs(y, get_data_from_file('projection_0008_map_16.xvg'), 'black', '--')
	line1, = PDFs(y,get_data_from_file('3130_DSPC_projection.xvg'),'CG DSPC', '#d62d20', '--')
	handles.append(line1)
	plt.xlabel('PC projection value (A)')
	plt.ylabel('Probability density (a.u.)')
	
	plt.legend(handles = handles)
	plt.show()

	print("K-S-S coeff: %s, %s" % (stats.ks_2samp((AA_f), (MAP_f))))

if __name__ == "__main__":
	main()
