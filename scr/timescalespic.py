import matplotlib.pyplot as plt


def get_data_from_file(file_name, time):
	data = []	
	file_ = open(file_name, 'r')
	lines = file_.readlines()
	file_.close()
	for line in lines:
		if line.find('E**2') != -1:
			data_1 = lines[lines.index(line) + 1:lines.index(line)+ 11]
		if line.find('E**1')!=-1:
			data_2 = lines[lines.index(line) + 1:lines.index(line)+ 11]
	if time == 't2':
		return [float(line.split()[1]) for line in data_1] 
	elif time == 't1':
		return [float(line.split()[1]) for line in data_2]

def PDFs(data, label, col, style=''):
	p = plt.loglog([i + 1 for i in range(10)], data, label = label, color = col, linestyle = style)
	return p


def main(file1,file2,type_,time,fout):	
	handles = []
	if time == 't2':
		data = get_data_from_file(file1, time)
		line1, = PDFs(data, r'First trajectory, $\tau_2$', 'blue', '--')
		handles.append(line1)
		data = get_data_from_file(file2, time)
		line1, = PDFs(data, r'Second trajectory, $\tau_2$','green', '--')
		handles.append(line1)
	elif time == 't1':
		data = get_data_from_file(file1, time)
		line1, = PDFs(data, r'First trajectory, $\tau_1$', 'blue', '--')
		handles.append(line1)
		data = get_data_from_file(file2, time)
		line1, = PDFs(data, r'Second trajectory, $\tau_1$','green', '--')
		handles.append(line1)
	plt.ylim([.1, 4 * 10 ** 2])
	plt.xlim([1, 10])
	plt.ylabel('Relaxation time (ns)')
	plt.xlabel('Component')
	plt.legend(handles = handles, ncol = 2)
	if type_ == 'auto':
		plt.title('Autocorrelation relaxation')
		plt.savefig(fout)
		print('Picture saved as '+fout)
	elif type_ == 'kss':
		plt.title('KSS relaxation')
		plt.savefig(fout)
		print('Picture saved as '+fout)
	
