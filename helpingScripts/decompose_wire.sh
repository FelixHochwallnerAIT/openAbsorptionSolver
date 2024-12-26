#!/bin/bash

# Input
case_folder=./workingFolder/wire

# Switch to case folder
cd $case_folder

# Remove processor folders if they exist
rm -rf processor*

# Get time directories from periodicFlow directory
mv periodicFlow/* .

# Get all time directories
timeDirs=$(find . -maxdepth 1 -type d -name '[0-9]*' -printf "%f\n")

# Sort time directories
timeDirs=$(echo $timeDirs | tr " " "\n" | sort -n)

# Save first and last value of timeDirs to a new variable
firstTimeDir=$(echo $timeDirs | awk '{print $1}')
lastTimeDir=$(echo $timeDirs | awk '{print $NF}')

# Decompose case
# decomposePar -allRegions -time  > log.decomposePar 2>&1
decomposePar -allRegions -time $firstTimeDir:$lastTimeDir \
    > log.decomposePar 2>&1

# Remove 0 time directory
timeDirs=$(echo $timeDirs | sed 's/0//')

# Get processor directories
processorDirs=$(find . -maxdepth 1 -type d -name 'processor*')

# Loop over processor directories
# and create periodic flow directory
# and then copy all time directories into this directory
for processorDir in $processorDirs; do
    mkdir $processorDir/periodicFlow
    for timeDir in $timeDirs; do
        mv $processorDir/$timeDir $processorDir/periodicFlow/
    done
done

# Move back all timeDirs to periodicFlow directory
for timeDir in $timeDirs; do
    mv $timeDir periodicFlow/
done