#!/usr/bin/python
import glob, os, shutil, random
from subprocess import call

def writeToFasta(filename, seq_name, seq):
	print(filename)
	file = open(filename, "w")
	file.write(seq_name + "\n")
	
	db_file = filename
	spec_file = filename.split(".")[0] + ".msalign"


	with open("local_test_two_ptm.sh", "a") as sh:
		sh.write("toppic -d -k " + db_file + " " + spec_file + "\n")
	
	for i in range(0, len(seq)):
		file.write(seq[i])
		if (i % 70 == 69):
			file.write("\n")
		

def processOneFile(filename):

	lines =  open(filename + ".fasta").read().splitlines()

	print("two_ptm/" + filename[5:])

	if not os.path.exists("two_ptm/" + filename[5:]): 
		os.mkdir("two_ptm/" + filename[5:])
		
	seq_name = lines[0]
	seq = ""
	for i in range(1, len(lines)):
		seq = seq + lines[i]
	
	i = 0
	j = 0
	
	seq_bak = seq
	k = 0
	
	for i in range(1, len(seq) - 21):
		seq = seq_bak
		if seq[i] == "T":
			for j in range(i + 20, len(seq) - 1):
				seq = seq_bak
				if seq[j] == "S":
					k = k + 1
					if k <= 3:
						seq = seq[0:i] + "B" + seq[i + 1:j] + "X" + seq[j + 1:]
						fname = "two_ptm/" + filename[5:] + "/" + filename[5:] + "-" + str(i) + "-" + str(j) + ".fasta"
						writeToFasta(fname, seq_name, seq)
						shutil.copyfile(filename + ".msalign", "two_ptm/" + filename[5:] + "/" + filename[5:] + "-" + str(i) + "-" + str(j) + ".msalign")
				


def main():
	
	if not os.path.exists("two_ptm"): 
		os.mkdir("two_ptm")
			
	with open("local_test_two_ptm.sh", "a") as sh:
		sh.write("#!/bin/bash\n")
	
	files = glob.glob("data/*.fasta")
	for f in files:
		processOneFile(f.split(".")[0])
	

if __name__ == "__main__":
   main()

