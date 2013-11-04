#!/usr/bin/env python
'''
This script is used to benchmark the program.
Square lattice with deifferent parameters is used.
'''
from __future__ import division
__author__ = "Vladimir Iglovikov"

import os
import time

#create simlink of the ggeom program in the 
try:
	os.symlink(os.path.join(os.getcwd(), "..", "..", "EXAMPLE", "geom", "ggeom"), os.path.join(os.getcwd(), "ggeom"))
except:
	pass

inList = ["small_td0.in", "small_td1.in", "median_td0.in", "median_td1.in", "large_td0.in", "large_td1.in"]

fName = open("results", "w")
fName.close()

resultTimes = {}
for inFile in inList:
	startTime = time.time()
	os.system("./ggeom {inFile} > {logFile}".format(inFile=inFile, logFile = "log" + inFile ))
	resultTimes[inFile] = time.time() - startTime
	fName = open("results", "a")
	print >> fName, inFile, " used ", resultTimes[inFile], "seconds"	
	fName.close()
	print inFile, " used ", resultTimes[inFile], "seconds"	
