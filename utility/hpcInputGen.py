#!/usr/bin/python

import os
import sys
import shutil
import time
import subprocess
import glob
import numpy as np
import math

from itertools import imap
from operator import sub
from operator import add
from os import mkdir
from os import chdir
from os import getcwd
from os import remove
from os import listdir
from shutil import copy
from numpy import array 

# Usage: nebGen.py initial final NumberOfImages

def readAAcoord(fname):
    with open(fname,'r') as inpos:
        data = inpos.readlines()
    atomType       = []
    atomNumber     = 0
    atomTypeNumber = 0
    atomAsi        = []
    atomConstraint = []
    atomicNumbers  = []
    atomComment    = []
    latticeInfo    = []
    shift          = 0

    try:
        atomTypeNumber  = int(data[0].split()[0])
    except:
        print 'WARNING: IGNORE THIS UNLESS LATTICE CONSTANT DATA IS INCLUDED IN FORT.970 and FORT.971'
        shift = 5
        latticeInfo=data[0:5]
        atomTypeNumber  = int(data[shift].split()[0])

    try:
        atomNumber = int(data[shift].split()[1])
    except:
        sys.exit('ERROR: The number of atoms must be integer.')

    try:
        atomicNumbers = [int(x) for x in data[shift].split()[2:2+atomTypeNumber]]
    except:
        sys.exit('ERROR: The number of atomic numbers is not given or it is not integer.')

    for line in data[1+shift:1+shift+atomNumber]:
	items = line.split()
        atomType.append(items[0])
	atomAsi.append(np.array([float(x) for x in items[1:4]]))

        if len(items) == 4:
            sys.exit('ERROR: please write atomic constraint type in the 5th column')

        if len(items) > 4  and items[4] != '/':
            atomConstraint.append(items[4])
            atomComment.append('/')            

    return latticeInfo,atomType,atomTypeNumber,atomNumber,atomAsi,atomConstraint,atomicNumbers,atomComment


def writeFort970(fname,latticeInfo,atomTypeNumber,atomType,atomNumber,atomAsi,atomConstraint,atomicNumbers,atomComment):
    outpos = open(fname,'w')
    
    if not latticeInfo:
        pass
    else:
        for i in range(len(latticeInfo)):
            outpos.write(latticeInfo[i])

    outpos.write("{0:<3d} {1:4<d} ".format(atomTypeNumber,atomNumber))
    [outpos.write("{0:>4d} ".format(int(x))) for x in atomicNumbers ]
    outpos.write("\n")

    for i in range(len(atomAsi)):
	outpos.write("{0:<3s} {1:<20.16f} {2:<20.16f} {3:<20.16f} {4:<3s}".format(atomType[i], atomAsi[i][0], atomAsi[i][1], atomAsi[i][2], atomConstraint[i]))
	[outpos.write(" {0:<s} ".format(comment)) for comment in atomComment[i]]
	outpos.write("\n")

def writeVnorm(fname,atomAsi):
    outpos = open(fname,'w')
    
    for i in range(len(atomAsi)):
	outpos.write("{0:<20.16f} {1:<20.16f} {2:<20.16f}".format(atomAsi[i][0], atomAsi[i][1], atomAsi[i][2]))
	outpos.write("\n")

class Collective:
    def __init__(self):
        self.latticeInfo     = []
        self.atomTypeNumber  = []
        self.atomType        = []
        self.atomNumber      = []
        self.atomAsi         = []
        self.atomConstraint  = []
        self.atomicNumbers   = []
        self.atomComment     = []
    def appending(self,tmpLatticeInfo,tmpAtomType,tmpAtomTypeNumber,tmpAtomNumber,tmpAtomAsi,tmpAtomConstraint,tmpAtomicNumbers,tmpAtomComment):
        self.latticeInfo.append(tmpLatticeInfo)
        self.atomType.append(tmpAtomType)
        self.atomTypeNumber.append(tmpAtomTypeNumber)
        self.atomNumber.append(tmpAtomNumber)
        self.atomAsi.append(tmpAtomAsi)
        self.atomConstraint.append(tmpAtomConstraint)
        self.atomicNumbers.append(tmpAtomicNumbers)
        self.atomComment.append(tmpAtomComment)

def checkerDiff(classCollective):
    if classCollective.atomType[0] != classCollective.atomType[1]:
        print 'ERROR:atomType'
    if (classCollective.atomNumber[0] != classCollective.atomNumber[1]).any():
        print 'ERROR:atomNumber'
        print classCollective.atomNumber[0]
        print classCollective.atomNumber[1]
    if (classCollective.atomicNumbers[0] != classCollective.atomicNumbers[1]).any():
        print 'ERROR:atomicNumbers'
        print classCollective.atomicNumbers[0]
        print classCollective.atomicNumbers[1]
    if classCollective.atomConstraint[0] != classCollective.atomConstraint[1]:
        print 'ERROR:atomConstraint'
    if classCollective.atomComment[0] != classCollective.atomComment[1]:
        print 'ERROR:atomComment'
    tmp =  [np.subtract(x2,x1) for x2, x1 in zip(collective.atomAsi[1],collective.atomAsi[0])]
    i = 0
    for posarray in tmp:
        if posarray[0] > 0.8:        # means initial shifted to left on x dir.
            print 'working00', i, posarray[0]
            print 'initial & final = ', collective.atomAsi[0][i][0], collective.atomAsi[1][i][0]
            collective.atomAsi[0][i][0] += 1.0
            print 'initial & final = ', collective.atomAsi[0][i][0], collective.atomAsi[1][i][0]
        elif posarray[0] < -0.8:
            print 'working01', i, posarray[0]
            print 'initial & final = ', collective.atomAsi[0][i][0], collective.atomAsi[1][i][0]
            collective.atomAsi[1][i][0] += 1.0
            print 'initial & final = ', collective.atomAsi[0][i][0], collective.atomAsi[1][i][0]
        if posarray[1] > 0.8:
            print 'working10', i, posarray[1]
            print 'initial & final = ', collective.atomAsi[0][i][1], collective.atomAsi[1][i][1]
            collective.atomAsi[0][i][1] += 1.0
            print 'initial & final = ', collective.atomAsi[0][i][1], collective.atomAsi[1][i][1]
        elif posarray[1] < -0.8:
            print 'working10', i, posarray[1]
            print 'initial & final = ', collective.atomAsi[0][i][1], collective.atomAsi[1][i][1]
            collective.atomAsi[1][i][1] += 1.0
            print 'initial & final = ', collective.atomAsi[0][i][1], collective.atomAsi[1][i][1]
        if posarray[2] > 0.8:
            print 'working20', i, posarray[2]
            print 'initial & final = ', collective.atomAsi[0][i][0], collective.atomAsi[1][i][2]
            collective.atomAsi[0][i][2] += 1.0
            print 'initial & final = ', collective.atomAsi[0][i][0], collective.atomAsi[1][i][2]
        elif posarray[2] < -0.8:
            print 'working21', i, posarray[2]
            print 'initial & final = ', collective.atomAsi[0][i][0], collective.atomAsi[1][i][2]
            collective.atomAsi[1][i][2] += 1.0
            print 'initial & final = ', collective.atomAsi[0][i][0], collective.atomAsi[1][i][2]
        i += 1

outputName = []
initial   = sys.argv[1]
final     = sys.argv[2]

try:
    numImage = int(sys.argv[3])
except:
    sys.exit('ERROR: Number of the images must be integer')

tmpAtomNumber     = 0
tmpAtomTypeNumber = 0
tmpAtomConstraintNumber = 0
tmpAtomAsi        = []
tmpAtomConstraint = []
tmpAtomComment    = []

# ------ reading and putting into list inside class -------

collective = Collective()

tmpLatticeInfo,tmpAtomType,tmpAtomTypeNumber,tmpAtomNumber,tmpAtomAsi,tmpAtomConstraint,tmpAtomicNumbers,tmpAtomComment = readAAcoord(initial)
collective.appending(tmpLatticeInfo,tmpAtomType,tmpAtomTypeNumber,tmpAtomNumber,tmpAtomAsi,tmpAtomConstraint,tmpAtomicNumbers,tmpAtomComment)

tmpLatticeInfo,tmpAtomType,tmpAtomTypeNumber,tmpAtomNumber,tmpAtomAsi,tmpAtomConstraint,tmpAtomicNumbers,tmpAtomComment = readAAcoord(final)
collective.appending(tmpLatticeInfo,tmpAtomType,tmpAtomTypeNumber,tmpAtomNumber,tmpAtomAsi,tmpAtomConstraint,tmpAtomicNumbers,tmpAtomComment)

# ----- error checking -----
#checkerDiff(collective)

# ----- printing -----

if numImage == 0:
    sys.exit('The number of the images must be more than 0')
elif numImage != 0:
    tmp0 = map(sub,collective.atomAsi[1],collective.atomAsi[0])
    for i in range(1,numImage+1):
        tmp1 = [x*i/(numImage+1) for x in tmp0]
        tmp2 = map(add,collective.atomAsi[0],tmp1)
        collective.atomAsi.append(tmp2)

collective.atomAsi.append(collective.atomAsi.pop(1))

k = 0
for atomAsiSet in collective.atomAsi:
    i = 0
    for atomAsiData in atomAsiSet:
        if atomAsiData[0] >= 1.0:
            print 'working 0', k, i
            atomAsiData[0] -= 1.0
        if atomAsiData[1] >= 1.0:
            print 'working 1', k, i
            atomAsiData[1] -= 1.0
        if atomAsiData[2] >= 1.0:
            print 'working 2', k, i
            atomAsiData[2] -= 1.0
        i += 1
    k += 1

for i in range(2 + numImage):
    if i <= 9:
        mkdir('0'+str(i))
        chdir('0'+str(i))
    elif i > 9:
        mkdir(str(i))
        chdir(str(i))
    writeVnorm('fort.973',tmp0)
    outputName.append('fort.970')
    writeFort970(outputName[i],collective.latticeInfo[0],collective.atomTypeNumber[0],collective.atomType[0],collective.atomNumber[0],collective.atomAsi[i],collective.atomConstraint[0],collective.atomicNumbers[0],collective.atomComment[0])
    chdir('../')
