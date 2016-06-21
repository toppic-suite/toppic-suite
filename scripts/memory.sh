#!/bin/sh
while [ true ]
do
  ps -p $1 -o %cpu,%mem
  sleep 20
done
