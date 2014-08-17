#!/usr/bin/python

import sys, getopt

def getLine(line, clines):
	arr = line.split("\t")
	num = int(arr[0])
	res = ""
	for i in range(1, len(clines)):
		cline = clines[i]
		arr = cline.split("\t")
		
		if (int(arr[4]) == num):
			res = cline
			break
	
	arr = res.split("\t")
	return arr[4] + "\t" + arr[5] + "\t" + arr[6] + "\t" + arr[7] + "\t" + arr[9] + "\t" + arr[11] + "\t" + arr[17] + "\t" + arr[18] + "\n"


def checkIn(line, clines):
	res = 0
	arr1 = line.split("\t")
	for i in range(1, len(clines)):
		cline = clines[i]
		arr2 = cline.split("\t")
		if ((arr1[3] == arr2[4]) and (arr1[9] == arr2[11]) ):
			if float(arr1[15]) > float(arr2[17]):
				print arr1[3] + "\t" + arr1[9] + "\t" + arr1[15] + "\t" + arr2[17] 
				res = 1
			else: 
				res = 2
			break
	return res

def process(cppResult, javaResult, outfile):
	file = open(outfile, "w")

	with open(cppResult, 'r+') as c:
		clines = c.readlines()
		
	with open(javaResult, 'r+') as j:
		jlines = j.readlines()
		
	file.write("status" + "\t" + jlines[0].rstrip('\n') + "\t")
	arr = clines[0].split("\t")
	file.write(arr[4] + "\t" + arr[5] + "\t" + arr[6] + "\t" + arr[7] + "\t" + arr[9] + "\t" + arr[11] + "\t" + arr[17] + "\t" + arr[18] + "\n")
	for i in range(1, len(jlines)):
		line = jlines[i]
		if (int(checkIn(line, clines)) == 0):
			file.write("Java" + "\t" + line )
		#elif (int(checkIn(line, clines)) == 1):
		#	file.write("same" + "\t" + line.rstrip('\n') + "\t" + getLine(line, clines))
		elif (int(checkIn(line, clines)) == 1):
			file.write("Not same good protein" + "\t" + line )
			
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
