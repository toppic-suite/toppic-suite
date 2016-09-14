#!/bin/bash
## usage: bash copyright.sh ./base/

SCRIPT_PATH=$(pwd)

cd $1
for file in *.cpp; do
    cat $SCRIPT_PATH/copyright_template.txt > tempfile; 
    echo $file
    cat $file >> tempfile;
    mv tempfile $file
done

for file in *.hpp; do
    cat $SCRIPT_PATH/copyright_template.txt > tempfile; 
    echo $file
    cat $file >> tempfile;
    mv tempfile $file
done
