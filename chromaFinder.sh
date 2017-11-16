#!/bin/bash
# # File Preprocessing script
# Input
#	1. Bam file to be processed
#	2. Reference file
#	3. analysis files on or off
#	4. Only itration 1 (1 if only iteration i is needed 0 otherwise)
#	5. chromaclique threshold
#	6. subsection on or off
#	7. Chromosome Id as present in the sam file
#	8. chromosome Start
#	9. chromosome end 



cd ./

if [ $6 = "0" ]; then
  echo "Subsection Turned Off"
  echo "Process is going to take quite some time since the whole alignment file has to be checked. Please be patient."
  #3.Cleaning the bam file
  echo "Starting bamcleaner"
  bamCleaner $1 cleaned.bam
  echo "End bamcleaner"
else 
  echo "Subsection Turned On"
  #2.Cutting out a subsection from the supplied bam file
  echo "Starting samtools subsection"
  samtools view -b -h -o subsection.bam $1 $7:$8-$9 
  #Sorting the subsection of the Bam File
  #samtools sort -n -o subSorted.bam subsection.bam
  echo "End samtools subsection"

  #3.Cleaning the subsection bam file
  echo "Starting bamcleaner"
  bamCleaner subsection.bam cleaned.bam
  echo "End bamcleaner"
fi


#4.Sorting the cleaned bam file by name
echo "Starting samtools sort by name"
samtools sort -n -o cleanedSorted.bam cleaned.bam
echo "End samtools sort ny name"


#5.Checking the sorted bam file
echo "Starting bamchecker"
bamChecker cleanedSorted.bam ChromacliqueWorking.bam
echo "End bamcleaner"

#Removing all other by product files
rm subsection.bam
echo "Removed subsection.bam"
rm cleanedSorted.bam
echo "cleanedSorted.bam"

echo "File Preparation finished successfully"

if [ $4 = "0" ]; then
  echo "Starting chromaclique"
  chromaclique --analysisFiles --bam --uniqueSupportFilter=-1 --nome --no_singletons --nomeParam=$5 --minOverlap=2 --FS=./FSSwitches.txt --RS=./RSSwitches.txt --ref=$2 ./ChromacliqueWorking.bam ./Result
else
  chromaclique --iterations=1 --analysisFiles --bam --uniqueSupportFilter=-1 --nome --no_singletons --nomeParam=$5 --minOverlap=2 --FS=./FSSwitches.txt --RS=./RSSwitches.txt --ref=$2 ./ChromacliqueWorking.bam ./Result
fi

rm cleaned.bam
echo "Removed cleaned.bam"

if [ $3 = "0" ]; then
  rm ChromacliqueWorking.bam
  rm ResultReads.txt
  echo "Removed ChromacliqueWorking.bam"
fi

echo "Chromaclique ended successfully"












