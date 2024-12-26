#!/bin/bash

# Input
case_folder=./workingFolder/wire

# Switch to case folder
cd $case_folder

# Post process the steady state solution

# Get all time directories
timeDirs=$(find . -maxdepth 1 -type d -name '[0-9]*' -printf "%f\n")

# Sort time directories
timeDirs=$(echo $timeDirs | tr " " "\n" | sort -n)

# Get latest time directory
latestTimeDir=$(echo $timeDirs | awk '{print $NF}')

# Get all regions
regions=($(find constant -maxdepth 1 -type d -exec basename {} \;))
# Remove leading constant from regions
regions=("${regions[@]/constant/}")

# Filter out empty and non-directory elements
filtered_regions=()
for region in "${regions[@]}"; do
    if [[ -n $region && -d "constant/$region" ]]; then
        filtered_regions+=("$region")
    fi
done

# Backup 0 folder
mkdir -p backup
cp -r 0 backup
cp -r $latestTimeDir backup
mv postProcessing/ backup/postProcessing_steadyState

# Move steady state solution to 0 folder
for region in "${filtered_regions[@]}"; do
    mv $latestTimeDir/$region/* 0/$region/
done

# Remove all time directories but 0
rm -rf [1-9]* 0.[0-9]*

# Change fvSchemes ddtScheme steadyState to Euler in all regions
for region in "${filtered_regions[@]}"; do
    sed -i "s/\<steadyState\>/Euler/" "system/$region/fvSchemes"
done

# Remove relaxationFactors cw and T to 0.5 in all regions
for region in "${regions[@]}"; do
    sed -i '/fields/,/}/{//!d}' system/$region/fvSolution
    sed -i '/fields/a\  {' system/$region/fvSolution
done

# Change back controlDict
mv system/controlDict.bak system/controlDict
