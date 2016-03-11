#!/bin/python
import sys
import os
from decimal import Decimal

print "msalign file: " + os.path.abspath(sys.argv[1])

msalign_file = open(os.path.abspath(sys.argv[1]), 'r').readlines()

count = 0

i = 1

while i < len(msalign_file):
	tmp = msalign_file[i].split("\t")
	if int(tmp[17]) >= 10:
		count = count + 1
	i = i + 1

print count
