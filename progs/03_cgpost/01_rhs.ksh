#!/bin/bash
# Read in the first x betas, 
# The correct prefix must be supplied as an argument
# for example, for Maryland use mdcg
#
if [ "$#" -ne 1 ]
then
     echo "rhs: Incorrect number of arguments."
     echo "Usage: rhs prefix"
     echo "Usage: prefix is everything before .coef"
     exit 1
fi

#Store the prefix
prefix=$1
#
# Get the number of betas
numrec=`awk 'BEGIN {  RS = " " } {print $1  }' ${prefix}.coef | tail +4 | wc -l` 
echo $numrec
#
# Create the data file
head -n $numrec ${prefix}.betas | awk ' BEGIN { ORS = " " } {print $2}' > 01_rhs.dat
#
# Get the length of the record
length=`cat 01_rhs.dat | wc -c`
#
# Create the SAS file
echo "libname temp '.';" > 01_rhs.sas
echo " " >> 01_rhs.sas
echo "data temp.rhs;" >> 01_rhs.sas
echo "infile '01_rhs.dat' lrecl=$length;" >>01_rhs.sas;
echo "input" >> 01_rhs.sas
awk 'BEGIN {  RS = " " } {print $1  }' ${prefix}.coef | tail +4 >> 01_rhs.sas
echo ";" >> 01_rhs.sas
echo " " >> 01_rhs.sas
echo "proc contents;" >> 01_rhs.sas
echo "run;" >> 01_rhs.sas
echo " " >> 01_rhs.sas
echo "proc print;" >> 01_rhs.sas
echo "run;" >> 01_rhs.sas
#
# Run the program
sas 01_rhs.sas
# Clean Up
#rm -f 01_rhs.dat
#rm -f 01_rhs.sas
#rm -f 01_rhs.log
