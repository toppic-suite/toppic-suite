#!/usr/bin/python
import glob, os, shutil
from subprocess import call

def rmChar(s):
	s = s[s.find('.')+1:]
	s = s[0:s.find('.')]
	return s

def writeToFasta(filename, seq_name, seq):
	file = open(filename, "w")
	file.write(seq_name + "\n")
	
	for i in range(0, len(seq)):
		file.write(seq[i])
		if (i % 70 == 69):
			file.write("\n")
		

def processOneLine(line):
	if not os.path.exists("data"): 
		os.mkdir("data")
	lines = line.split("\t")
	filename = "data/spec" + lines[4] + ".fasta"
	seq_name = lines[11]
	seq = lines[15]
	seq = rmChar(seq)
	writeToFasta(filename, seq_name, seq)
	shutil.move("spec" + lines[4] + ".msalign", "data/spec" + lines[4] + ".msalign")
	
def main():
	flag = False
	files = open("ec_hcd.msalign").read().splitlines()
	#for f in files:
	for i in range(0, len(files)):
		f = files[i]
		if (f.startswith("SCANS=")):
			flag = True
			scan_num = f[6:]
			file_name = "spec" + f[6:] + ".msalign"
			with open(file_name, "a") as spec:
				spec.write("BEGIN IONS\n")
				spec.write(files[i - 1] +"\n")
				spec.write(files[i] +"\n")
		elif (f.startswith("END")):
			with open(file_name, "a") as spec:
				spec.write(f + "\n")
				spec.close()
			flag = False
		elif (flag):
			with open(file_name, "a") as spec:
				spec.write(f + "\n")			
	
	lines =  open("ec_hcd_no_ptm").read().splitlines()
	scan = []
	for i in range(1,len(lines)):
		processOneLine(lines[i])
		scan.append(int(lines[i].split("\t")[4]))
	##print(scan)



if __name__ == "__main__":
   main()
