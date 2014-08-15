#!/usr/bin/python

import sys, getopt


def writeToFile(file, line):
	file.write(line)        
	
	
def process(inputfile):
	file_ptm = open(inputfile + ".PTM", "w")
	with open(inputfile) as fp:
		for line in fp:
			arr = line.split("\t")
			if line.startswith("Data_file_name"):
				writeToFile(file_ptm, line)
				continue
				
			if int(arr[14]) > 0:
				writeToFile(file_ptm, line)
					
	file_ptm.close()

def main(argv):
	if len(sys.argv) == 1 :
		print 'find_ptm_ids.py -i <inputfile> '
		sys.exit()
	try:
		opts, args = getopt.getopt(argv,"hi:",["ifile="])
	except getopt.GetoptError:
		print 'find_ptm_ids.py -i <inputfile> '
		sys.exit()
	for opt, arg in opts:
		if opt == '-h':
			print 'find_ptm_ids.py -i <inputfile> '
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		else:
			print 'find_ptm_ids.py -i <inputfile> '
   
	print 'Input file is:\t', inputfile
	
	process(inputfile)

if __name__ == "__main__":
   main(sys.argv[1:])
