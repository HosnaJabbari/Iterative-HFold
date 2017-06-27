#!/bin/csh

if ($#argv != 2) then
    echo "Usage: $0 <software> <directory>"
    echo "   <software> can be simfold"
    echo "      Creating separate packages for pairfold and multifold is NOT supported any more!"
    echo "      Use the whole MultiRNAFold package if you are interested in pairfold or multifold."
    echo "Example:"
    echo "$0 simfold ../simfold"
    #echo "$0 pairfold ../pairfold"
    #echo "$0 multifold ../multifold" 
    exit    
endif    

set soft = $argv[1]
set dir = $argv[2]


#if (!(($soft =~ "simfold") ||($soft =~ "pairfold") || ($soft =~ "multifold"))) then
if (!($soft =~ "simfold")) then
    echo "Usage: $0 <software> <directory>"
    echo "   <software> can be simfold"
    echo "      Creating separate packages for pairfold and multifold is NOT supported any more!"
    echo "      Use the whole MultiRNAFold package if you are interested in pairfold or multifold."
    echo "Example:"
    echo "$0 simfold ../simfold"
    exit
endif    

if (-e "$dir") then
    echo "Cannot create dir: $dir already exists."
    exit
endif     
 
echo "Making directory $dir ..."
mkdir $dir

echo "Copying parameters ..."
cp -rf params/ $dir

echo "Copying sources ..."
mkdir $dir/src
cp -rf src/common/ $dir/src/
cp -rf src/$soft/ $dir/src/

echo "Copying header files ..."
mkdir $dir/include
cp -rf include/init.h include/constants.h include/externs.h include/globals.h include/structs.h include/$soft.h $dir/include

echo "Copying driver ..."
cp -rf $soft*.cpp $dir
if ($soft == "simfold") then
    cp test_get_counts.cpp $dir
endif


echo "Copying Makefile ..."
cp -rf Makefile_$soft $dir/Makefile

echo "Done."

