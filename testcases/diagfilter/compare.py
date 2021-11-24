#!/usr/bin/python

import sys, getopt

def checkIn(line, clines):
	res = 0
	arr1 = line.split("\t")
	for i in range(1, len(clines)):
		cline = clines[i]
		arr2 = cline.split("\t")
		if ((arr1[3] == arr2[1]) and (arr1[9] == arr2[6].rstrip('\n'))):
			res = 1
			break
	return [res, i]

def process(cppResult, javaResult, outfile):
	file = open(outfile, "w")

	with open(cppResult, 'r+') as c:
		clines = c.readlines()
		
	with open(javaResult, 'r+') as j:
		jlines = j.readlines()
		
	file.write("status" + "\t" + jlines[0].rstrip('\n') + "\t")
	file.write(clines[0])
	for i in range(1, len(jlines)):
		line = jlines[i]
		status = checkIn(line, clines)
		if (status[0] == 0):
			file.write("Java" + "\t" + line.rstrip('\n') + "\t" + clines[status[1]])
		#elif (status[0] == 1):
		#	file.write("same" + "\t" + line.rstrip('\n') + "\t" + clines[status[1]])
			
	file.close()

	
def main(argv):
	if len(sys.argv) == 1 :
		print 'compare.py -c <cpp> -j <java> -o <output>'
		sys.exit()
	try:
		opts, args = getopt.getopt(argv,"hc:j:o:",["cpp=","java=", "output="])
	except getopt.GetoptError:
		print 'compare.py -c <cpp> -j <java> -o <output>'
		sys.exit()
	for opt, arg in opts:
		if opt == '-h':
			print 'compare.py -c <cpp> -j <java> -o <output>'
			sys.exit()
		elif opt in ("-c", "--cpp"):
			cppResult = arg
		elif opt in ("-j", "--java"):
			javaResult = arg
		elif opt in ("-o", "--output"):
			outfile = arg
		else:
			print 'compare.py -c <cpp> -j <java> -o <output>'
   
	print 'C++ result:\t', cppResult
	print 'Java result:\t', javaResult
	print 'Result:\t', outfile
	process(cppResult, javaResult, outfile)
	

if __name__ == "__main__":
   main(sys.argv[1:])
