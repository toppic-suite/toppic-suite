#!/usr/bin/python

import sys, getopt

def checkIn(line, clines):
	res = 0
	arr1 = line.split("\t")
	for i in range(1, len(clines)):
		cline = clines[i]
		arr2 = cline.split("\t")
		if ((arr1[0] == arr2[4]) and (arr1[4] == arr2[11]) ):
			res = 1
			break
		elif (arr1[0] == arr2[0]):
			res = 2
			break
	return res

def process(cppResult, javaResult):
	file = open("res.txt", "w")

	with open(cppResult, 'r+') as c:
		clines = c.readlines()
		
	with open(javaResult, 'r+') as j:
		jlines = j.readlines()
		
	file.write("status" + "\t" + clines[0])
	
	for i in range(1, len(jlines)):
		line = jlines[i]
		if (int(checkIn(line, clines)) == 0):
			file.write("Java" + "\t" + line )
		elif (int(checkIn(line, clines)) == 1):
			file.write("same" + "\t" + line )
		elif (int(checkIn(line, clines)) == 2):
			file.write("Not same protein" + "\t" + line )
			
	file.close()

	
def main(argv):
	if len(sys.argv) == 1 :
		print 'check.py -c <cpp> -j <java>'
		sys.exit()
	try:
		opts, args = getopt.getopt(argv,"hc:j:",["cpp=","java="])
	except getopt.GetoptError:
		print 'check.py -c <cpp> -j <java>'
		sys.exit()
	for opt, arg in opts:
		if opt == '-h':
			print 'check.py -c <cpp> -j <java>'
			sys.exit()
		elif opt in ("-c", "--cpp"):
			cppResult = arg
		elif opt in ("-j", "--java"):
			javaResult = arg
		else:
			print 'check.py -c <cpp> -j <java>'
   
	print 'C++ result:\t', cppResult
	print 'Java result:\t', javaResult
	process(cppResult, javaResult)
	

if __name__ == "__main__":
   main(sys.argv[1:])
