#!/bin/python

from Bio import SeqIO

import itertools

fasta_file = "human_h3.fasta"

N = 120

NPTM = 3

AA = "ARNDCEQGHILKMFPSTWYVUO"

ptm_map = {}
ptm_map['ph'] = {'S', 'T', 'Y'}
ptm_map['acetyl'] = {'K', 'R'}
#ptm_map['meth'] = {'K', 'R'}
#ptm_map['dimeth'] = {'K', 'R'}
#ptm_map['trimeth'] = {'R'}

pos_list = list()

for p in ptm_map.values():
	pos_list = pos_list + list(p)

aa_list = list(set(pos_list))

aa_map = {}

for a in aa_list:
	aa_map[a] = pos_list.count(a)

def get_num(seq):
	all_sum = 0
	for s in list(itertools.combinations(seq, NPTM)):
		tmp = 1
		for p in s:
			tmp = tmp * aa_map[p]
		all_sum = all_sum + tmp
	return all_sum

def gen_all_pos(seq):
	seq_str = ''
	for s in seq:
		if s in aa_list:
			seq_str = seq_str + s
	return get_num(seq_str)

def gen_all_substr(seq, n):
	all_sum = 0
	for i in range(0, len(seq) - n + 1):
		all_sum = all_sum + gen_all_pos(seq[i:(i+n)])

	return all_sum

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

for fasta in fasta_sequences:
	print gen_all_substr(fasta.seq.tostring(), N)
