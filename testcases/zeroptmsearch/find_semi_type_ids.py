#!/usr/bin/python

import sys, getopt

def getStatus(pep):
	if (pep.startswith("M.") or pep.startswith(".")):
		if (pep.endswith(".")):
			return 1
		else:
			return 2
	else:
		if (pep.endswith(".")):
			return 3
		else:
			return 4

def writeToFile(file, arr):
	file.write(arr[3] + "\t" + arr[4] + "\t" + arr[5] + "\t")
	file.write(arr[8] + "\t" + arr[9] + "\t" + arr[13] + "\t")
	file.write(arr[15] + "\t" + arr[16] + "\n")
	
	
def process(inputfile, num):
	file1 = open(inputfile + ".COMPLETE", "w")
	file2 = open(inputfile + ".PREFIX", "w")
	file3 = open(inputfile + ".SUFFIX", "w")
	file4 = open(inputfile + ".INTERNAL", "w")
	with open(inputfile) as fp:
		for line in fp:
			arr = line.split("\t")
			if line.startswith("Data_file_name"):
				writeToFile(file1, arr)
				writeToFile(file2, arr)
				writeToFile(file3, arr)
				writeToFile(file4, arr)
				continue
				
			if int(arr[14]) == num:
				status = getStatus(arr[13])
				if status == 1:
					writeToFile(file1, arr)
				elif status == 2:
					writeToFile(file2, arr)
				elif status == 3:
					writeToFile(file3, arr)
				elif status == 4:
					writeToFile(file4, arr)
					
	file1.close()
	file2.close()
	file3.close()
	file4.close()

def main(argv):
	if len(sys.argv) == 1 :
		print 'find_semi_type_ids.py -i <inputfile> -n <PTM number>'
		sys.exit()
	try:
		opts, args = getopt.getopt(argv,"hi:n:",["ifile=","num="])
	except getopt.GetoptError:
		print 'find_semi_type_ids.py -i <inputfile> -n <PTM number>'
		sys.exit()
	for opt, arg in opts:
		if opt == '-h':
			print 'find_semi_type_ids.py -i <inputfile> -n <PTM number>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-n", "--num"):
			num = int(arg)
		else:
			print 'find_semi_type_ids.py -i <inputfile> -n <PTM number>'
   
	print 'Input file is:\t', inputfile
	print 'PTM num is:\t', num
	
	process(inputfile, num)

if __name__ == "__main__":
   main(sys.argv[1:])
