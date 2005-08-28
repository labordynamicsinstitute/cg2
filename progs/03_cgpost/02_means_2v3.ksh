#!/bin/bash
# Read in the first x means, 
# The correct prefix must be supplied as an argument
# for example, for Maryland use mdcg
#
if [ "$#" -ne 1 ]
then
     echo "means: Incorrect number of arguments."
     echo "Usage: means prefix"
     echo "Usage: prefix is everything before .coef"
     exit 1
fi

#Store the prefix
prefix=$1
#
# Get the number of means
numrec=`awk 'BEGIN {  RS = " " } {print $1  }' ${prefix}.coef | tail +4 | wc -l` 
numrec=$(($numrec+1))
echo $numrec
#
# Create the data file
head -n $numrec ${prefix}.means | awk ' BEGIN { ORS = " " } {print $1}' > 02_means_2v3.dat
#
# Get the length of the record
length=`cat 02_means_2v3.dat | wc -c`
#
# Create the SAS file
echo "libname temp '.';" > 02_means_2v3.sas
echo " " >> 02_means_2v3.sas
echo "data temp.means02;" >> 02_means_2v3.sas
echo "infile '02_means_2v3.dat' lrecl=$length;" >>02_means_2v3.sas;
echo "input" >> 02_means_2v3.sas
awk 'BEGIN {  RS = " " } {print $1  }' ${prefix}.coef | head -1 >> 02_means_2v3.sas
awk 'BEGIN {  RS = " " } {print $1  }' ${prefix}.coef | tail +4 >> 02_means_2v3.sas
echo ";" >> 02_means_2v3.sas
echo " " >> 02_means_2v3.sas
echo "proc contents;" >> 02_means_2v3.sas
echo "run;" >> 02_means_2v3.sas
echo " " >> 02_means_2v3.sas
echo "proc print;" >> 02_means_2v3.sas
echo "run;" >> 02_means_2v3.sas
#
# Run the program
sas 02_means_2v3.sas
# Clean Up
#rm -f 02_means_2v3.dat
#rm -f 02_means_2v3.sas
#rm -f 02_means_2v3.log
