#!/usr/bin/python
import glob, os, shutil, random
from subprocess import call

def writeToFasta(filename, seq_name, seq):
	print(filename)
	file = open(filename, "w")
	file.write(seq_name + "\n")
	
	db_file = filename
	spec_file = filename.split(".")[0] + ".msalign"


	with open("local_test_one_ptm.sh", "a") as sh:
		sh.write("toppic -d -k " + db_file + " " + spec_file + "\n")
	
	for i in range(0, len(seq)):
		file.write(seq[i])
		if (i % 70 == 69):
			file.write("\n")
		

def processOneFile(filename):

	lines =  open(filename + ".fasta").read().splitlines()

	print("one_ptm/" + filename[5:])

	if not os.path.exists("one_ptm/" + filename[5:]): 
		os.mkdir("one_ptm/" + filename[5:])
		
	seq_name = lines[0]
	seq = ""
	for i in range(1, len(lines)):
		seq = seq + lines[i]
	
	i = 0
	j = 0
	
	seq_bak = seq
	
	for i in range(1, len(seq) - 1):
		seq = seq_bak
		if seq[i] in "STY":
			if seq[i] == 'S':
				seq = seq[0:i] + "X" + seq[i + 1:]
			if seq[i] == 'T':	
				seq = seq[0:i] + "B" + seq[i + 1:]
					
			fname = "one_ptm/" + filename[5:] + "/" + filename[5:] + "-" + str(i) + ".fasta"
			writeToFasta(fname, seq_name, seq)
			shutil.copyfile(filename + ".msalign", "one_ptm/" + filename[5:] + "/" + filename[5:] + "-" + str(i) + ".msalign")
				


def main():
	
	if not os.path.exists("one_ptm"): 
		os.mkdir("one_ptm")
			
	with open("local_test_one_ptm.sh", "a") as sh:
		sh.write("#!/bin/bash\n")
	
	files = glob.glob("data/*.fasta")
	for f in files:
		processOneFile(f.split(".")[0])
	

if __name__ == "__main__":
   main()

