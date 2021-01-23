#!/usr/bin/env python

import csv
import os
import re
import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument('--n', help='Enter the number of vertices')
parser.add_argument('--outputFile', help='Enter the output seeds file path')

args = parser.parse_args()
n = int(args.n)
outputFilePath = args.outputFile
print "Creating partition 1 vertices : {}".format(n/2)

outputFile = open(outputFilePath, 'w')
i = 0
while i < n:
    val = random.randint(1, 1001)
    if val%2 == 1:
        outputFile.write(str(i)+"\n")
    i += 1

outputFile.close()
print "Partion created at {}".format(outputFilePath)
