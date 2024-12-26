#!/bin/bash

# Input
case_folder=./workingFolder/wire

# Switch to case folder
cd $case_folder

# Simulate the case steady state

# Get all time directories
timeDirs=$(find periodicFlow -maxdepth 1 -type d -name '[0-9]*' -printf "%f\n")

# Sort time directories
timeDirs=$(echo $timeDirs | tr " " "\n" | sort -n)

# Get first time directory
firstTimeDir=$(echo $timeDirs | awk '{print $1}')

# Copy first time directory
cp -r periodicFlow/$firstTimeDir .

# Backup controlDict
cp system/controlDict system/controlDict.bak

# Change controlDict
# Change startTime to firstTimeDir
sed -i "s/startTime .*;/startTime $firstTimeDir;/" system/controlDict
# Change endTime to 1000
sed -i "s/endTime .*;/endTime 1000;/" system/controlDict
# Change deltaT to 1
sed -i "s/deltaT.*;/deltaT 1;/" system/controlDict
# Change writeControl to timeStep
sed -i "0,/writeControl/s/writeControl.*;/writeControl timeStep;/" system/controlDict
# Change writeInterval to 100
sed -i "s/writeInterval.*;/writeInterval 100;/" system/controlDict

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

# Change fvSchemes ddtScheme Euler to steadyState in all regions
for region in "${filtered_regions[@]}"; do
    sed -i "s/\<Euler\>/steadyState/" "system/$region/fvSchemes"
done

# Change relaxationFactors cw and T to 0.5 in all regions
for region in "${regions[@]}"; do
    sed -i '/fields/,/}/{//!d}' system/$region/fvSolution
    sed -i '/fields/a\  {\n    T 0.5;\n    cw 0.5;' system/$region/fvSolution
done
