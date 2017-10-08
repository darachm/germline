#!/bin/bash

for i in $*
#For each file given as an argument to this script, set it to i
do 
  infile=$(pwd)'/'$i
#This variable is the complete path to the file
  outfile=$(pwd)'/'$( echo -n $i | sed 's/.dv/.png/' )
#This variable is the complete path for the new output, 
#and swaps .dv for .png
  fileBase=$(basename $infile)
#This variable strips the extension off
  echo $fileBase;
#Report what we're working on
#The below command calls imagej, and runs the expression in the
#quotes behind -e
  imagej -e '
run("Bio-Formats Importer", "open='$infile' color_mode=Default view=[Data Browser] stack_order=XYCZT use_virtual_stack");
run("Z Project...", "projection=[Max Intensity]");
run("Split Channels");
selectWindow("C1-MAX_'$fileBase'");
run("Enhance Contrast","saturated=0.3")
selectWindow("C2-MAX_'$fileBase'");
run("Enhance Contrast","saturated=0.3")
run("Merge Channels...", "c1=C2-MAX_'$fileBase' c3=C1-MAX_'$fileBase' create");
run("Stack to RGB");
saveAs("PNG", "'$outfile'");
run("Quit");
  '
done;
