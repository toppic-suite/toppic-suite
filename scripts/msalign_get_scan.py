#!/bin/python
import sys
import os
from decimal import Decimal
from util2 import *

if len(sys.argv) < 4:
    print "Usage: " + sys.argv[0] + " msalign_file_name scan"
    sys.exit()

print "msalign file: " + os.path.abspath(sys.argv[1])

scans = [int(i) for i in sys.argv[3:]]

print "Scans: " + str(scans)

msalign_file = open(os.path.abspath(sys.argv[1]), 'r').readlines()

spec_lst = readMsAlign(os.path.abspath(sys.argv[1]))

print "In total: " + str(len(spec_lst)) + " spectra"

for i in xrange(len(spec_lst)):
    if spec_lst[i].scan in scans:
		spec_lst[i].to_file(sys.argv[2])

