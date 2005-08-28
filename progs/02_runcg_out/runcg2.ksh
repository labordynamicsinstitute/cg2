#!/bin/bash
# Set up CG
#
# CG RUN Vars
CG=/home/kmckinney/cg2f90/bin/cg2_4v3
#CG=/data/working/02/creec003/fortpgms/cgf90/readbin/cg2/cg2_4v3
CGPARMFILE=cg.cgin
CGLOGFILE=cg.cglog
RUNCGLOG=runcg.log
BINFILE=cgout4
BETAFILE=cg.betas
MEANFILE=cg.means
STUB=cg
PROGDIR=/home/kmckinney/runcg8

# Clean up previous runs
rm -f ${STUB}.cgin
rm -f ${STUB}.cglog
rm -f ${STUB}.bin4
rm -f ${STUB}.betas
rm -f ${STUB}.means

# Get the number of variables
echo "%include '${PROGDIR}/config_param.sas';
         %bincalc(4);
" |
sas -stdio 2>setup.tmp
NUMVAR=$(grep "num coef" setup.tmp | awk '{ print $4 }')
((NUMVAR=NUMVAR-1))
/bin/rm -f setup.tmp

# Create the CGIN file
grep "cg.cgin" ${PROGDIR}/02_size_calc.log | tail -1 | awk '{ print $2 " " $3 " " $4 " " $5 " " var }' var=$NUMVAR > $CGPARMFILE

# No links needed for this version

rm -f ${STUB}.bin4
ln -s $BINFILE ${STUB}.bin4

    if [[ -e ${CGPARMFILE} ]]
    then
        echo "Running CG program ......................"
        echo "Starting CG program " >> ${RUNCGLOG}
        echo "(Log file of CG algorithm in ${CGLOGFILE}) " >> ${RUNCGLOG}
        date >> ${RUNCGLOG}
        $CG ${STUB} >> ${RUNCGLOG}
        echo "CG program has finished " >> ${RUNCGLOG}
        date >> ${RUNCGLOG}
        echo " " >> ${RUNCGLOG}
        echo "Regression summary:" >> ${RUNCGLOG}
        head -n 5 ${CGLOGFILE} >> ${RUNCGLOG}
        echo " " >> ${RUNCGLOG}
        grep "Total memory" ${CGLOGFILE} >> ${RUNCGLOG}
        tmp=`tail -n 1 ${CGLOGFILE}`
        set -- ${tmp}
        echo "Final Iteration         "  "${2}" >> ${RUNCGLOG}
        shift 2
        echo "$*" >> ${RUNCGLOG}
    else
        echo "Cannot run CG because of missing parameter file."
            echo "${rc_failed}"
        echo "Cannot run CG because of missing parameter file." >> ${RUNCGLOG}
        echo "Exiting." >> ${RUNCGLOG}
        exit 2
    fi
