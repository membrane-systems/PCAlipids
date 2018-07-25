def main(filename):
	file = open(filename, 'r')
	proj_n = 1
	file_out = open('proj_1.xvg', 'w')
	for line in file:
		if line.find('@') != -1:
			continue
		elif line.find('&') != -1:
			file_out.close()
			proj_n += 1
			file_out = open('proj_%s.xvg' % proj_n, 'w')
		else:
			file_out.write(line.split()[1] + '\n')
	file_out.close()
	file.close()