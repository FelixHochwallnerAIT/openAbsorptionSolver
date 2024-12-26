#!/bin/bash

# Input
case_folder=./workingFolder/wire

# Switch to case folder
cd $case_folder

# Prepare case for heatMassTransferWireFoam

# ----------------------
# Move all time directories (but 0) to periodicFlow directory
# and create periodicFlow directory if it doesn't exist

# Get all time directories
timeDirs=$(find . -maxdepth 1 -type d -name '[0-9]*' -printf "%f\n")

# Sort time directories
timeDirs=$(echo $timeDirs | tr " " "\n" | sort -n)

# Remove 0 time directory
timeDirs=$(echo $timeDirs | sed 's/0//')

# Create periodicFlow directory if it doesn't exist
if [ ! -d "periodicFlow" ]; then
    mkdir periodicFlow
fi

# Move all time directories to periodicFlow directory
for timeDir in $timeDirs; do
    mv $timeDir periodicFlow
done

# Path to the periodicFlow folder
periodicFlow_folder="periodicFlow"

# Path to the periodicFlow OpenFOAM dict file
dict_file="constant/periodicFlow"

# Create the periodicFlow file if it doesn't exist
if [ ! -f ${dict_file} ]; then
    touch ${dict_file}
fi

# Get the foldernames in the periodicFlow folder
foldernames=($(ls -d ${periodicFlow_folder}/*/))

# Remove trailing slashes from foldernames
foldernames=("${foldernames[@]%/}")

# Echo the foldernames
# echo "foldernames: ${foldernames[@]}"

# Create the timesteps array
timesteps=()
for foldername in "${foldernames[@]}"; do
    timesteps+=("${foldername##*/}")
done

# Echo the timesteps array
# echo "timesteps: ${timesteps[@]}"

# Write header to OpenFOAM dict file
echo "/*--------------------------------*- C++ -*----------------------------------*\\" > ${dict_file}
echo "  =========                 |" >> ${dict_file}
echo "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox" >> ${dict_file}
echo "   \\\\    /   O peration     | Website:  https://openfoam.org" >> ${dict_file}
echo "    \\\\  /    A nd           | Version:  7" >> ${dict_file}
echo "     \\\\/     M anipulation  |" >> ${dict_file}
echo "\\*---------------------------------------------------------------------------*/" >> ${dict_file}

# Write timesteps array to OpenFOAM dict file
echo "FoamFile" >> ${dict_file}
echo "{" >> ${dict_file}
echo "    version     2.0;" >> ${dict_file}
echo "    format      ascii;" >> ${dict_file}
echo "    class       dictionary;" >> ${dict_file}
echo "    location    \"constant\";" >> ${dict_file}
echo "    object      periodicFlow;" >> ${dict_file}
echo "}" >> ${dict_file}
echo "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" >> ${dict_file}
echo "timesteps" >> ${dict_file}
echo "(" >> ${dict_file}

# Loop through timesteps array and write each timestep to the dict file
for timestep in "${timesteps[@]}"; do
    echo "    \"${timestep}\"" >> ${dict_file}
done

# Close the timesteps array in the dict file
echo ");" >> ${dict_file}
echo "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" >> ${dict_file}

# Echo completion message
echo "Timesteps written to ${dict_file}"
