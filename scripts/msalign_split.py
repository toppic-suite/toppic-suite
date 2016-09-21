#!/bin/python
import sys
import os
from decimal import Decimal
from util2 import *

if len(sys.argv) < 3:
    print "Usage: " + sys.argv[0] + " msalign_file_name split_number"
    sys.exit()

print "msalign file: " + os.path.abspath(sys.argv[1])

msalign_file = open(os.path.abspath(sys.argv[1]), 'r').readlines()
prefix = os.path.abspath(sys.argv[1])[0:(len(os.path.abspath(sys.argv[1])) - 8)]

spec_lst = spec_lst = readMsAlign(os.path.abspath(sys.argv[1]))

print "In total: " + str(len(spec_lst)) + " spectra"
print "split into " + sys.argv[2] + " files"
suffix = 0
t = int(round(len(spec_lst) / float(sys.argv[2]) + 0.5))
for i in xrange(len(spec_lst)):
    if i % t == 0:
        suffix = suffix + 1
    spec_lst[i].to_file(str(prefix) + "_" + str(suffix) + ".msalign")
