#Steps to Run the grouping program
#1 Run firmcells.sas to create the worker-firm cells sorted by SEIN PIK
#2 Create the links below

rm -f fort.4 fort.7 fort.8 groups.in groups.log

ln -s groups      fort.4
ln -s ../cellsout fort.7
ln -s firmcells   fort.8

#3 Create the groups program input parameters
#cat ../cg.cgin | awk '{ print $2 " " $3 " " $4 }' > groups.in
ln -s ../cg.cgin groups.in

# Run the grouping program
/home/kmckinney/cg2f90/bin/groups < groups.in > groups.log 
