#!/usr/bin/env python

import csv
import os
import re
import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument('--n', help='Enter the number of vertices')
parser.add_argument('--outputFile', help='Enter the output seeds file path')
parser.add_argument(
    '--seedsPercent', help='Enter the expected percent of seeds')

args = parser.parse_args()
n = int(args.n)
outputFilePath = args.outputFile
seedsPercent = float(args.seedsPercent)
chance = 30
numberOfSeeds = int(n/seedsPercent)

print "Number of Seeds : {}".format(numberOfSeeds)

currentNumberOfSeeds = 0
iter = 0
isSeed = [False]*n

while currentNumberOfSeeds < numberOfSeeds:
    i = 0
    while i < n:
        val = random.randint(1, chance)
        if val == chance-1:
            if isSeed[i] == False:
                currentNumberOfSeeds += 1
                isSeed[i] = True
                if currentNumberOfSeeds == numberOfSeeds:
                    break
        i += 1
    iter += 1


# print "Number of iterations : {}".format(iter)


outputFile = open(outputFilePath, 'w')
i = 0
while i < n:
    if isSeed[i] == True:
        outputFile.write(str(i)+"\n")
    i += 1

outputFile.close()
