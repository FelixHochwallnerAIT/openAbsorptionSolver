#!/bin/bash

# Directories
sourceDir=./workingFolder/wire_init
targetDir=./workingFolder/wire_Antoine

# Print the directories
echo "Source directory: $sourceDir"
echo "Target directory: $targetDir"
echo

# Get time directories of source case
sourceTimeDirs=$(find $sourceDir \
    -maxdepth 1 \
    -type d \
    -name "[0-9]*" \
    -printf "%f\n")

# Sort the time directories
sourceTimeDirs=$(echo $sourceTimeDirs | tr " " "\n" | sort -n)

# Remove the 0 time directory
sourceTimeDirs=$(echo $sourceTimeDirs | sed 's/0//')

# Check if mesh of source case is 2D, if so make it 3D
checkMeshSource=$(checkMesh -case "$sourceDir" 2>&1)
checkMeshTarget=$(checkMesh -case "$targetDir" -region air 2>&1)
sourceDimensions=$(echo "$checkMeshSource" \
    | grep "geometric (non-empty/wedge) directions" \
    | awk '{print $3}')
sourceHeight=$(echo "$checkMeshSource" \
    | grep "bounding box" \
    | awk '{print $10}')
# Remove trailing ) from height
sourceHeight=$(echo $sourceHeight | sed 's/)//')
targetHeight=$(echo "$checkMeshTarget" \
    | grep "bounding box" \
    | awk '{print $10}')
# Remove trailing ) from height
targetHeight=$(echo $targetHeight | sed 's/)//')

if [ $sourceDimensions == 2 ]; then
    echo "Source case is 2D. Transforming to 3D..."
    # Change empyt boundary in blockMeshDict
    sed -i "s/empty/patch/g" $sourceDir/system/blockMeshDict
    # Recreate mesh
    blockMesh -case $sourceDir >> $sourceDir/log.blockMesh 2>&1
    echo "blockMesh changed."
    # Change empty boundary in all time folders
    for timeDir in $sourceTimeDirs; do
        sed -i "s/empty/zeroGradient/g" $sourceDir/$timeDir/p
        sed -i "s/empty/zeroGradient/g" $sourceDir/$timeDir/U
        rm "$sourceDir/$timeDir/total(p)"
        rm "$sourceDir/$timeDir/phi"
    done
    echo "Boundaries in time folders for p and U changed."
    factor=$(echo $targetHeight / $sourceHeight | bc -l)
    transformPoints -case $sourceDir -scale "(1 1 $factor)" \
        >> $sourceDir/log.transformPoints 2>&1
    echo "Mesh extruded."
    # Check if 2D mesh transformation was successful
    checkMeshSource=$(checkMesh -case "$sourceDir" 2>&1)
    sourceDimensions=$(echo "$checkMeshSource" \
        | grep "geometric (non-empty/wedge) directions" \
        | awk '{print $3}')
    sourceHeight=$(echo "$checkMeshSource" \
        | grep "bounding box" \
        | awk '{print $10}')
    # Remove trailing ) from height
    sourceHeight=$(echo $sourceHeight | sed 's/)//')
    targetHeight=$(echo "$checkMeshTarget" \
        | grep "bounding box" \
        | awk '{print $10}')
    # Remove trailing ) from height
    targetHeight=$(echo $targetHeight | sed 's/)//')
    if [ $targetHeight != $sourceHeight ]; then
        echo "ERROR: Meshes have different heights now. Aborting."
        exit 1
    fi
    if [ $sourceDimensions != 3 ]; then
        echo "ERROR: Mesh transformation failed. Aborting."
        exit 1
    fi
    echo "Mesh transformation successful!"
elif [ $sourceDimensions == 3 ]; then
    echo "Source case is 3D."
    # Check if meshes have the same height 
    if [ $targetHeight != $sourceHeight ]; then
        echo "Meshes have different heights. Aborting."
        exit 1
    fi
else
    echo "Source case is neither 2D nor 3D. Aborting."
    exit 1
fi
echo "Meshes have the same height."

# Get regions of target cast
REGIONS=$(find $targetDir/constant \
    -mindepth 1 \
    -maxdepth 1 \
    -type d \
    -printf "%f\n")

# Remove the cell coordinates file for all regions in target case
for region in $REGIONS; do
    rm -f $targetDir/0/$region/C*
    rm -f $targetDir/0/$region/cellDist
done

# Remove phi from air
rm -f $targetDir/0/air/phi

for timeDir in $sourceTimeDirs; do
    echo "Mapping fields for time $timeDir"
    
    # Create the target time directory
    cp -r $targetDir/0 $targetDir/$timeDir

    # Change startTime
    sed -i "s/startTime.*;/startTime $timeDir;/g" \
        $targetDir/system/controlDict

    # Change startFrom to startTime
    sed -i "s/startFrom.*;/startFrom startTime;/g" \
        $targetDir/system/controlDict

    mapFields \
        -case $targetDir \
        -consistent \
        -sourceTime $timeDir \
        -targetRegion air \
        $sourceDir \
        >> $targetDir/log.mapFields_air 2>&1
done

echo "All fields have been mapped."

# Change back controlDict
# Set startTime to 0
sed -i "s/startTime.*;/startTime 0;/g" \
    $targetDir/system/controlDict

# Set startFrom to latestTime
sed -i "s/startFrom.*;/startFrom latestTime;/g" \
    $targetDir/system/controlDict