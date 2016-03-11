#!/usr/bin/python
import glob, os, shutil
from subprocess import call
import operator, shutil

def writeToFasta(filename, seq_name, seq):
	file = open(filename, "w")
	file.write(seq_name + "\n")
	
	for i in range(0, len(seq)):
		file.write(seq[i])
		if (i % 70 == 69):
			file.write("\n")

### S->Z
def processFasta(fname):
	#print fname.split(".")[0]
	lines = open("data/seq/" + fname).read().splitlines()
	seq_name = lines[0]
	seq = reduce(operator.add, lines[1:])
	count = 0
	sh_file = open("test1.sh", "a")
	##print seq
	for i in xrange(len(seq)):
		if seq[i] == 'S':
			count = count + 1
			new_seq = seq[0:i] + 'Z' + seq[i + 1:]
			##print new_seq
			shutil.copyfile("data/" + fname.split(".")[0] + ".msalign", "data/seq1/" + fname.split(".")[0] + ".msalign" + str(count))
			writeToFasta("data/seq1/" + fname + str(count), seq_name, new_seq)
			sh_file.write("./test_graph -i testcases/mass_graph/mass_graph_mods.txt -k -p 0 -j 20 ")
			sh_file.write("data/seq1/" + fname + str(count) + " ")
			sh_file.write("data/seq1/" + fname.split(".")[0] + ".msalign" + str(count) + "\n")


if not os.path.exists("data"): 
	os.mkdir("data")

if not os.path.exists("data/seq1"): 
	os.mkdir("data/seq1")

fasta_lst = os.listdir("data/seq")

for f in fasta_lst:
	processFasta(f)
