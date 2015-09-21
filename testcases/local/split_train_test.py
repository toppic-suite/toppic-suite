#!/usr/bin/python
import sys, random 

def split(fname):
    files = open(fname).read().splitlines()
    n = len(files) - 1
    random.seed(42)
    train_lst = random.sample(xrange(n), n / 2)
    train_lst.sort()
    f_train = open(fname + ".train", "w")
    f_test = open(fname + ".test", "w")
    f_train.write(files[0] + "\n")
    f_test.write(files[0] + "\n")
    for i in range(0, n):
        if (i in train_lst):
            f_train.write(files[1 + i] + "\n")
        else:
            f_test.write(files[1 + i] + "\n")

if __name__ == "__main__":
	split(sys.argv[1])

