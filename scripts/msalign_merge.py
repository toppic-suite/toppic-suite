#!/bin/python
import sys
import os
from decimal import Decimal
from util2 import *

N = 100000

if len(sys.argv) < 3:
    print "Usage: " + sys.argv[0] + " merged_msalign_file_name file_list"
    sys.exit()

os.remove(sys.argv[1])

file_list = [i for i in sys.argv[2:]]

print "File list: " + str(file_list)

for i in xrange(len(file_list)):
	print "Processing: ", os.path.abspath(file_list[i])
	spec_lst = readMsAlign(os.path.abspath(file_list[i]))
	for sp in spec_lst:
		sp.id = sp.id + i * N
		sp.scan = sp.scan + i * N
		sp.ms_one_id = sp.ms_one_id + i * N
		sp.ms_one_scan = sp.ms_one_scan + i * N
		sp.feature_id = sp.feature_id + i * N
		sp.to_file(sys.argv[1])

