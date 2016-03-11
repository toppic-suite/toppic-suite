#!/usr/bin/python
import sys
import os
import shutil
from decimal import Decimal

if len(sys.argv) < 3:
    print "Usage: " + sys.argv[0] + " file_name copy_number"
    sys.exit()

abs_path = os.path.abspath(sys.argv[1]) 
print "file name: " + abs_path

file_base_name, file_ext = os.path.splitext(abs_path)
print "base name: " + file_base_name

for i in range(1, int(sys.argv[2])):
    dest = file_base_name + "_" + str(i) +  file_ext
    print dest
    shutil.copyfile(abs_path, dest)
