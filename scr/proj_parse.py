def main(filename):
	file = open(filename, 'r')
	proj_n = 1
	name = filename.split('.')[0]
	file_out = open('%s_1.xvg' % name, 'w')
	for line in file:
		if line.find('@') != -1:
			continue
		elif line.find('&') != -1:
			file_out.close()
			proj_n += 1
			file_out = open('%s_%s.xvg' % (name, proj_n), 'w')
		else:
			file_out.write(line.split()[1] + '\n')
	file_out.close()
	file.close()