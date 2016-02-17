#!/bin/env python
# ccs = check crab status
import os

for root, dirs, files in os.walk("crab_projects"):
    for dirname in dirs:
	if(not(dirname == "results" or dirname == "inputs")): 
	    print "crab status -d crab_projects/"+dirname+" | tee -a /tmp/theofil/status.txt;"


print ""
print ""
print ""

for root, dirs, files in os.walk("crab_projects"):
    for dirname in dirs:
	if(not(dirname == "results" or dirname == "inputs")): 
	    print "crab report -d crab_projects/"+dirname+" | tee -a /tmp/theofil/reports.txt;"

print ""
print ""
print ""

for root, dirs, files in os.walk("crab_projects"):
    for dirname in dirs:
	if(not(dirname == "results" or dirname == "inputs")): 
	    print "cp crab_projects/"+dirname+"/results/lumiSummary.json /tmp/theofil/"+dirname+".json"
