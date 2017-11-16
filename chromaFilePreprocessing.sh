#!/bin/bash
# # File Preprocessing script
# Input
#	1. Bam file to be processed
#	2. Output file path
#	3. Reference file

#Generating the switch rate profiles
echo "Starting GCSwitchFinder"
GCSwitchfinder $1 $2 $3
echo "End bamcleaner" 
